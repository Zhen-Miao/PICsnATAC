"""Theoretical distributions and insertion rate estimation for PIC counts.

Translates R/Insertion_rates_and_DAR_test.R (estimation functions) to Python.
"""

import math

import numpy as np
from joblib import Parallel, delayed
from scipy.optimize import minimize_scalar
from scipy.special import comb

from picsnATAC._utils import parse_peak_names


def _get_theoretical_c1(insertion_rate, peak_length=500.0):
    """Compute theoretical PIC distribution under condition 1 (no size filtering).

    Single-allele distribution is computed from closed-form expressions,
    then convolved with itself for bi-allelic scenario.

    Parameters
    ----------
    insertion_rate : float
        Insertion rate per 1000 bp.
    peak_length : float
        Width of peak in bp.

    Returns
    -------
    np.ndarray
        Array of length 6: P(count=0) through P(count>=5), after bi-allelic
        convolution.
    """
    lam = insertion_rate * peak_length / 1000.0
    p = np.zeros(6)
    p[0] = np.exp(-0.5 * lam) * (2 - np.exp(-0.5 * lam))
    p[1] = np.exp(-0.5 * lam) * (lam - 2 + 2 * np.exp(-0.5 * lam))
    p[2] = np.exp(-lam) * (np.exp(0.5 * lam) * (lam**2 - 4 * lam + 8) - 8) / 4
    p[3] = np.exp(-lam) * (
        np.exp(0.5 * lam) * (lam**3 - 6 * lam**2 + 24 * lam - 48) + 48
    ) / 24
    p[4] = np.exp(-lam) * (
        np.exp(0.5 * lam) * (lam**4 - 8 * lam**3 + 48 * lam**2 - 192 * lam + 384)
        - 384
    ) / 192
    p[5] = 1 - np.sum(p[:5])

    # Bi-allelic convolution
    # R: convolve(p, rev(p), type="open") = np.convolve(p, p, mode='full')
    p_conv = np.convolve(p, p)
    p_conv[5] = np.sum(p_conv[5:11])
    return p_conv[:6]


def _insertion_to_c1(insertion_rates=None, peak_lengths=None):
    """Build insertion_rate × peak_length lookup matrix of expected counts.

    Parameters
    ----------
    insertion_rates : np.ndarray, optional
        Default: np.arange(0.01, 20.01, 0.01)
    peak_lengths : np.ndarray, optional
        Default: np.arange(200, 1050, 50)

    Returns
    -------
    np.ndarray
        Shape (n_rates, n_lengths).
    """
    if insertion_rates is None:
        insertion_rates = np.arange(1, 2001) * 0.01
    if peak_lengths is None:
        peak_lengths = np.arange(4, 21) * 50

    mat = np.zeros((len(insertion_rates), len(peak_lengths)))
    weights = np.arange(7, dtype=np.float64)  # 0,1,2,3,4,5,6
    for i, rate in enumerate(insertion_rates):
        for j, plen in enumerate(peak_lengths):
            p = _get_theoretical_c1(rate, plen)
            # p has 6 elements (indices 0-5), weights has 7 (0-6)
            # R: sum(p_W_s_theo * 0:6) where p_W_s_theo has length 6
            # This means indices 0..5 weighted by 0..5, and the last element
            # (p[5] = P(>=5)) is weighted by 5
            mat[i, j] = np.sum(p * weights[:6])
    return mat


def _get_theoretical_c12(insertion_rate, peak_length=1000.0, pad_length=0.0,
                         dirs=(-1, 1), count_type="PIC",
                         min_frag_length=25.0, max_frag_length=600.0,
                         cap_insertion=20):
    """Compute theoretical PIC distribution under conditions 1+2 (with size filtering).

    Parameters
    ----------
    insertion_rate : float
        Insertion rate per 1000 bp.
    peak_length : float
        Width of peak in bp.
    min_frag_length, max_frag_length : float
        Fragment size filtering bounds.
    cap_insertion : int
        Maximum number of insertions in one peak region.

    Returns
    -------
    np.ndarray
        Probability of observing 1, 2, ..., cap_fragment PIC counts
        (after bi-allelic convolution). Length = cap_insertion // 2.
    """
    lam = insertion_rate * peak_length / 1000.0
    cap_fragment = cap_insertion // 2
    p_W_m = np.zeros(cap_fragment)

    for mi in range(cap_fragment):
        m = mi + 1  # 1-indexed
        all_index_n = np.arange(m + 1, cap_insertion + 1)
        p_W_m_n = np.zeros(len(all_index_n))
        for ni, n in enumerate(all_index_n):
            val = (comb(n - 1, m, exact=True) * 0.5 ** (n - 1) *
                   lam ** n * np.exp(-lam)) / math.factorial(n)
            p_W_m_n[ni] = val
        p_W_m[mi] = np.sum(p_W_m_n)

    # Expected capture probability
    exp_capture = (
        np.exp(-min_frag_length * insertion_rate / 1000.0)
        - np.exp(-max_frag_length * insertion_rate / 1000.0)
    ) / (1 - np.exp(-peak_length * insertion_rate / 1000.0))

    p_W_s = np.zeros(cap_fragment)
    for ki in range(cap_fragment):
        k = ki + 1
        all_index_m = np.arange(k, cap_fragment + 1)
        p_W_s_k = np.zeros(len(all_index_m))
        for mii, m in enumerate(all_index_m):
            p_W_s_k[mii] = (comb(m, k, exact=True) *
                             exp_capture ** k *
                             (1 - exp_capture) ** (m - k) *
                             p_W_m[m - 1])
        p_W_s[ki] = np.sum(p_W_s_k)

    # Bi-allelic convolution
    p_W_s_full = np.concatenate([[1 - np.sum(p_W_s)], p_W_s])
    p_conv = np.convolve(p_W_s_full, p_W_s_full)
    # Return P(count=1) through P(count=cap_fragment)
    return p_conv[1:cap_fragment + 1]


def _insertion_to_c12(insertion_rates=None, peak_lengths=None,
                      pad_length=0.0, dirs=(-1, 1),
                      count_type="fragment",
                      min_frag_length=25.0, max_frag_length=600.0,
                      cap_insertion=20):
    """Build lookup matrix for condition 1+2 expected counts.

    Returns
    -------
    np.ndarray
        Shape (n_rates, n_lengths). Each entry is the expected count.
    """
    if insertion_rates is None:
        insertion_rates = np.arange(1, 2001) * 0.01
    if peak_lengths is None:
        peak_lengths = np.arange(4, 21) * 50

    mat = np.zeros((len(insertion_rates), len(peak_lengths)))
    for i, rate in enumerate(insertion_rates):
        for j, plen in enumerate(peak_lengths):
            p = _get_theoretical_c12(
                rate, plen, pad_length, dirs, count_type,
                min_frag_length, max_frag_length, cap_insertion,
            )
            # Expected count = sum(p[k] * (k+1))
            mat[i, j] = np.sum(p * np.arange(1, len(p) + 1))
    return mat


def _log_loss_frag_c1(est_inser, peak_length, capturing_rates, obs_pic_vec,
                      cap_fragment=5):
    """Log-likelihood under condition 1 with discrete capturing rates.

    Parameters
    ----------
    est_inser : float
        Estimated insertion rate per 1000 bp.
    peak_length : float
        Peak width.
    capturing_rates : np.ndarray
        Discretized capturing rate per cell.
    obs_pic_vec : np.ndarray
        Observed PIC counts per cell.
    cap_fragment : int
        Maximum count value.

    Returns
    -------
    float
        Log-likelihood (to be maximized).
    """
    p_theo = _get_theoretical_c1(est_inser, peak_length)

    all_rates = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99])
    ll = 0.0

    for ct in all_rates:
        # Build transition matrix: pt[observed, true_count]
        pt = np.zeros((cap_fragment + 1, cap_fragment + 1))
        for ot in range(cap_fragment + 1):
            for k in range(ot, cap_fragment + 1):
                pt[ot, k] = (comb(k, ot, exact=True) *
                             ct ** ot * (1 - ct) ** (k - ot) *
                             p_theo[k])
        ptv = pt.sum(axis=1)  # marginal probability of each observed count

        # Select cells with this capturing rate
        sel = capturing_rates == ct
        if not np.any(sel):
            continue
        obs_sel = obs_pic_vec[sel].astype(int)
        obs_sel = np.clip(obs_sel, 0, cap_fragment)
        probs = np.log(np.maximum(ptv[obs_sel], 1e-300))
        probs = np.clip(probs, -15, None)
        ll += np.sum(probs)

    return ll


def _log_loss_frag_ME(est_inser, peak_length, capturing_rates, obs_pic_vec,
                      min_frag_length=25.0, max_frag_length=600.0,
                      cap_fragment=10):
    """Log-likelihood under conditions 1+2 with per-cell capturing rates.

    Parameters
    ----------
    est_inser : float
        Estimated insertion rate.
    peak_length : float
        Peak width.
    capturing_rates : np.ndarray
        Per-cell capturing rates.
    obs_pic_vec : np.ndarray
        Observed PIC counts.
    min_frag_length, max_frag_length : float
        Size filtering bounds.
    cap_fragment : int
        Maximum count.

    Returns
    -------
    float
        Log-likelihood.
    """
    p_theo = _get_theoretical_c12(
        est_inser, peak_length,
        min_frag_length=min_frag_length,
        max_frag_length=max_frag_length,
    )
    # Prepend P(0)
    p0 = 1 - np.sum(p_theo)
    p_full = np.concatenate([[p0], p_theo])

    n_cells = len(obs_pic_vec)
    ll = 0.0

    for i in range(n_cells):
        ot = int(obs_pic_vec[i])
        ct = capturing_rates[i]
        total = 0.0
        for k in range(ot, cap_fragment + 1):
            if k < len(p_full):
                total += (comb(k, ot, exact=True) *
                          ct ** ot * (1 - ct) ** (k - ot) *
                          p_full[k])
        val = np.log(max(total, 1e-300))
        ll += max(val, -15)

    return ll


def obs_to_insertion_MLE_obj(pic_mat, capturing_rates, plen, n_cores=1):
    """Maximize log-likelihood per peak, return objective values.

    Parameters
    ----------
    pic_mat : scipy.sparse matrix or np.ndarray
        Peak-by-cell PIC count matrix.
    capturing_rates : np.ndarray
        Capturing rate per cell.
    plen : np.ndarray
        Peak lengths.
    n_cores : int
        Number of parallel jobs.

    Returns
    -------
    np.ndarray
        Max log-likelihood per peak.
    """
    if hasattr(pic_mat, "toarray"):
        pic_dense = pic_mat.toarray()
    else:
        pic_dense = np.asarray(pic_mat)
    n_peaks = len(plen)

    def _optimize_one(pp):
        res = minimize_scalar(
            lambda x: -_log_loss_frag_c1(x, plen[pp], capturing_rates, pic_dense[pp, :]),
            bounds=(0.01, 20),
            method="bounded",
        )
        return -res.fun

    results = Parallel(n_jobs=n_cores)(
        delayed(_optimize_one)(pp) for pp in range(n_peaks)
    )
    return np.array(results)


def obs_to_insertion_MLE_lam(pic_mat, capturing_rates, plen, n_cores=1):
    """Maximize log-likelihood per peak, return optimal insertion rates.

    Parameters
    ----------
    pic_mat : scipy.sparse matrix or np.ndarray
        Peak-by-cell PIC count matrix.
    capturing_rates : np.ndarray
        Capturing rate per cell.
    plen : np.ndarray
        Peak lengths.
    n_cores : int
        Number of parallel jobs.

    Returns
    -------
    np.ndarray
        Optimal insertion rate (per 1000 bp) per peak.
    """
    if hasattr(pic_mat, "toarray"):
        pic_dense = pic_mat.toarray()
    else:
        pic_dense = np.asarray(pic_mat)
    n_peaks = len(plen)

    def _optimize_one(pp):
        res = minimize_scalar(
            lambda x: -_log_loss_frag_c1(x, plen[pp], capturing_rates, pic_dense[pp, :]),
            bounds=(0.01, 20),
            method="bounded",
        )
        return res.x

    results = Parallel(n_jobs=n_cores)(
        delayed(_optimize_one)(pp) for pp in range(n_peaks)
    )
    return np.array(results)


def obs_to_insertion_ME(pic_mat, capturing_rates, cell_type_labels,
                        cap_insertion=20, min_frag_length=25.0,
                        max_frag_length=600.0, insertion_rates=None,
                        peak_lengths=None):
    """Moment estimator for insertion rates.

    Parameters
    ----------
    pic_mat : scipy.sparse matrix or np.ndarray
        Peak-by-cell PIC count matrix. Row names should be 'chr:start-end'.
    capturing_rates : np.ndarray
        Per-cell capturing rates.
    cell_type_labels : np.ndarray
        Cell type label per cell.
    cap_insertion : int
        Maximum insertions.
    min_frag_length, max_frag_length : float
        Size filtering bounds.
    insertion_rates, peak_lengths : np.ndarray, optional
        Ranges for lookup table.

    Returns
    -------
    np.ndarray
        Shape (n_peaks, n_cell_types) of estimated insertion rates.
    """
    if insertion_rates is None:
        insertion_rates = np.arange(1, 2001) * 0.01
    if peak_lengths is None:
        peak_lengths = np.arange(4, 21) * 50

    # Build lookup table
    c12_mat = _insertion_to_c12(
        insertion_rates=insertion_rates,
        peak_lengths=peak_lengths,
        min_frag_length=min_frag_length,
        max_frag_length=max_frag_length,
    )

    if hasattr(pic_mat, "toarray"):
        pic_dense = pic_mat.toarray()
    else:
        pic_dense = np.asarray(pic_mat)

    # Parse peak names to get peak lengths
    peak_names = getattr(pic_mat, "_peak_names", None)
    if peak_names is not None:
        _, pk_starts, pk_ends = parse_peak_names(peak_names)
    else:
        # Assume row names are stored externally - try to infer from matrix shape
        raise ValueError(
            "Peak names not available. Please pass a matrix with _peak_names attribute "
            "or provide peak_names separately."
        )

    plen = pk_ends - pk_starts

    # Group peaks by length bins
    plen_gap = int(peak_lengths[1] - peak_lengths[0])
    plen_group = np.clip(
        np.ceil((plen - peak_lengths[0]) / plen_gap).astype(int) + 1,
        1, len(peak_lengths),
    )

    # Clamp capturing rates
    cap_rates = np.maximum(capturing_rates, 0.2)

    cell_type_labels = np.asarray(cell_type_labels)
    cell_types = list(dict.fromkeys(cell_type_labels))

    n_peaks = pic_dense.shape[0]
    est_inser_all = np.zeros((n_peaks, len(cell_types)))

    for ct_idx, ct in enumerate(cell_types):
        sel = cell_type_labels == ct
        pic_sub = pic_dense[:, sel]
        cap_sub = cap_rates[sel]
        n_cells_ct = sel.sum()

        # w_bar = mean count adjusted for capturing rate
        w_bar = (pic_sub / cap_sub[np.newaxis, :]).sum(axis=1) / n_cells_ct

        unique_groups = np.unique(plen_group)
        for pg in unique_groups:
            pg_mask = plen_group == pg
            c12_col = c12_mat[:, pg - 1]  # 1-indexed group to 0-indexed column
            w_bar_sub = w_bar[pg_mask]

            # Find closest insertion rate for each peak
            diffs = np.abs(w_bar_sub[:, np.newaxis] - c12_col[np.newaxis, :])
            best_idx = diffs.argmin(axis=1)
            est_inser_all[pg_mask, ct_idx] = insertion_rates[best_idx]

    return est_inser_all
