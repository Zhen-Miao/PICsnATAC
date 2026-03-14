"""Differential Accessible Region (DAR) testing via Likelihood Ratio Test.

Translates the DAR_by_LRT function from R/Insertion_rates_and_DAR_test.R.
"""

import warnings

import numpy as np
from scipy import sparse
from scipy.stats import chi2

from picsnATAC._utils import parse_peak_names
from picsnATAC.insertion_rates import (
    _log_loss_frag_ME,
    obs_to_insertion_ME,
    obs_to_insertion_MLE_obj,
)


def DAR_by_LRT(pic_mat, capturing_rates, cell_type_labels,
               n_cores=1, plen=None, min_frag_length=25.0,
               max_frag_length=600.0, estimation_approach="MLE",
               peak_names=None):
    """Likelihood Ratio Test for Differential Accessible Regions.

    Tests for differential accessibility between exactly two cell types.

    Parameters
    ----------
    pic_mat : scipy.sparse matrix or np.ndarray
        Peak-by-cell PIC count matrix.
    capturing_rates : np.ndarray
        Capturing rate per cell.
    cell_type_labels : np.ndarray
        Cell type label per cell (must have exactly 2 unique types).
    n_cores : int
        Number of parallel jobs for MLE.
    plen : np.ndarray, optional
        Peak lengths. If None, parsed from peak_names or row names.
    min_frag_length, max_frag_length : float
        Fragment size filtering bounds (used in ME approach).
    estimation_approach : str
        Either 'MLE' or 'ME'.
    peak_names : list of str, optional
        Peak names in 'chr:start-end' format. Required if plen is None.

    Returns
    -------
    np.ndarray
        P-values, one per peak.
    """
    capturing_rates = np.asarray(capturing_rates, dtype=np.float64)
    cell_type_labels = np.asarray(cell_type_labels)

    # Validate inputs
    if len(capturing_rates) != len(cell_type_labels):
        raise ValueError("Number of cells do not match among inputs")

    n_cells_mat = pic_mat.shape[1] if sparse.issparse(pic_mat) else pic_mat.shape[1]
    if n_cells_mat != len(capturing_rates):
        raise ValueError("Number of cells do not match among inputs")

    if estimation_approach not in ("MLE", "ME"):
        raise ValueError("Please choose estimation_approach from MLE or ME")

    n_pks = pic_mat.shape[0]

    if n_pks * n_cells_mat >= 2**31 - 1:
        raise ValueError(
            "pic_mat too large, please slice the matrix by peaks and run the test"
        )

    ct_uniq = list(dict.fromkeys(cell_type_labels))
    if len(ct_uniq) != 2:
        raise ValueError("cell_type_labels must have exactly 2 unique cell types")

    # Get peak lengths
    if plen is None:
        if peak_names is None:
            raise ValueError("Either plen or peak_names must be provided")
        _, pk_starts, pk_ends = parse_peak_names(peak_names)
        plen = pk_ends - pk_starts
        plen = np.ceil(plen / 100).astype(int) * 100

    # Cap the counts
    if sparse.issparse(pic_mat):
        pic_mat = pic_mat.copy()
        cts = pic_mat.data.copy()
        cts[cts >= 10] = 0
        cts[cts > 5] = 5
        pic_mat.data[:] = cts
        pic_mat.eliminate_zeros()
        pic_dense = pic_mat.toarray()
    else:
        pic_dense = np.array(pic_mat, dtype=np.float64)
        pic_dense[pic_dense >= 10] = 0
        pic_dense[pic_dense > 5] = 5

    # Discretize capturing rates
    capturing_rates = np.ceil(capturing_rates * 10) / 10
    capturing_rates = np.clip(capturing_rates, 0.2, None)
    capturing_rates[capturing_rates > 0.9] = 0.99

    # Split by cell type
    mask_1 = cell_type_labels == ct_uniq[0]
    mask_2 = cell_type_labels == ct_uniq[1]

    pic_mat_1 = pic_dense[:, mask_1]
    pic_mat_2 = pic_dense[:, mask_2]
    cap_rates_1 = capturing_rates[mask_1]
    cap_rates_2 = capturing_rates[mask_2]

    if estimation_approach == "MLE":
        # Likelihood under null (pooled)
        ll_all = obs_to_insertion_MLE_obj(
            pic_dense, capturing_rates, plen, n_cores=n_cores,
        )
        # Likelihood under alternative (separate per cell type)
        ll_1 = obs_to_insertion_MLE_obj(
            pic_mat_1, cap_rates_1, plen, n_cores=n_cores,
        )
        ll_2 = obs_to_insertion_MLE_obj(
            pic_mat_2, cap_rates_2, plen, n_cores=n_cores,
        )
        test_stat = 2 * (ll_1 + ll_2 - ll_all)
        p_val = chi2.sf(test_stat, df=1)

    else:  # ME approach
        # Null: single group
        # Need to attach peak_names to matrix for obs_to_insertion_ME
        pic_with_names = pic_dense.copy()
        pic_with_names = np.asarray(pic_with_names)

        # We need to pass peak names through the ME function
        # Create a wrapper matrix-like object
        class _MatWithNames:
            def __init__(self, data, names):
                self._data = data
                self._peak_names = names
                self.shape = data.shape
            def toarray(self):
                return self._data
            def __getattr__(self, name):
                return getattr(self._data, name)

        mat_named = _MatWithNames(pic_dense, peak_names)

        lamb_all = obs_to_insertion_ME(
            mat_named, capturing_rates,
            cell_type_labels=np.full(len(capturing_rates), "A"),
            min_frag_length=min_frag_length,
            max_frag_length=max_frag_length,
        )
        lamb_all = lamb_all[:, 0]

        lamb_full = obs_to_insertion_ME(
            mat_named, capturing_rates,
            cell_type_labels=cell_type_labels,
            min_frag_length=min_frag_length,
            max_frag_length=max_frag_length,
        )
        lamb_full_1 = lamb_full[:, 0]
        lamb_full_2 = lamb_full[:, 1]

        ll_null = np.zeros(n_pks)
        ll_full = np.zeros(n_pks)

        for pki in range(n_pks):
            pl = plen[pki]
            ll_null[pki] = _log_loss_frag_ME(
                lamb_all[pki], pl, capturing_rates, pic_dense[pki, :],
                min_frag_length=min_frag_length,
                max_frag_length=max_frag_length,
            )
            ll_full[pki] = (
                _log_loss_frag_ME(
                    lamb_full_1[pki], pl, cap_rates_1, pic_mat_1[pki, :],
                    min_frag_length=min_frag_length,
                    max_frag_length=max_frag_length,
                ) +
                _log_loss_frag_ME(
                    lamb_full_2[pki], pl, cap_rates_2, pic_mat_2[pki, :],
                    min_frag_length=min_frag_length,
                    max_frag_length=max_frag_length,
                )
            )

        test_stat = 2 * (ll_full - ll_null)
        p_val = chi2.sf(test_stat, df=1)

    return p_val
