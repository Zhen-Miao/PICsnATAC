"""EM algorithm for estimating peak open probability and cell capturing rate.

Translates R/param_estimate_joint_p_q.R to Python.
"""

import numpy as np
from scipy import sparse


def get_r_by_ct_mat_pq(cell_type_set, r_by_c, cell_type_labels,
                        n_features_per_cell, p_acc=0.0005, q_acc=0.0005,
                        n_max_iter=800, verbose=True):
    """EM algorithm to estimate peak open probability and cell capturing rate.

    Parameters
    ----------
    cell_type_set : list of str
        Unique cell type labels.
    r_by_c : scipy.sparse matrix or np.ndarray
        Peak-by-cell matrix (will be binarized internally).
    cell_type_labels : np.ndarray
        Cell type label per cell.
    n_features_per_cell : int
        Number of peaks (rows).
    p_acc : float
        Convergence threshold for p.
    q_acc : float
        Convergence threshold for q.
    n_max_iter : int
        Maximum EM iterations.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        'p_by_t': np.ndarray (n_features, n_cell_types) - open probabilities
        'q_vec': np.ndarray (n_cells,) - capturing rates per cell
    """
    cell_type_labels = np.asarray(cell_type_labels)

    if sparse.issparse(r_by_c):
        col_names = getattr(r_by_c, '_col_names', None)
    else:
        col_names = None

    n_cells_total = r_by_c.shape[1]

    # Binarize: set all nonzero values to 1
    if sparse.issparse(r_by_c):
        r_by_c = r_by_c.copy().astype(np.float64)
        r_by_c.data[:] = 1.0
    else:
        r_by_c = (np.asarray(r_by_c) > 0).astype(np.float64)

    # Output containers
    p_by_t = np.zeros((n_features_per_cell, len(cell_type_set)))
    q_vec = np.zeros(n_cells_total)

    for gg_idx, gg in enumerate(cell_type_set):
        sel = cell_type_labels == gg
        if sparse.issparse(r_by_c):
            r_sub = r_by_c[:, sel]
            n_reads_in_cell = np.asarray(r_sub.sum(axis=0)).ravel()
            n_reads_in_region = np.asarray(r_sub.sum(axis=1)).ravel()
        else:
            r_sub = r_by_c[:, sel]
            n_reads_in_cell = r_sub.sum(axis=0)
            n_reads_in_region = r_sub.sum(axis=1)

        n_cell_sub = sel.sum()

        # Initialize
        q_cur = n_reads_in_cell / max(n_reads_in_cell.max(), 1)
        p_cur = n_reads_in_region / n_cell_sub

        diff1 = 1.0
        diff2 = 1.0
        n_iter = 0

        while (diff1 > p_acc or diff2 > q_acc) and n_iter < n_max_iter:
            p_old = p_cur.copy()
            q_old = q_cur.copy()

            # E-step: estimate q from p
            q_cur = n_reads_in_cell / max(p_cur.sum(), 1e-10)
            q_cur = np.minimum(q_cur, 0.999)

            # M-step: estimate p from q
            p_cur = n_reads_in_region / max(q_cur.sum(), 1e-10)
            p_cur = np.minimum(p_cur, 0.999)

            n_iter += 1

            # Convergence check
            with np.errstate(divide="ignore", invalid="ignore"):
                diff1 = np.nansum(np.abs(p_cur - p_old) / np.abs(p_cur)) / n_features_per_cell
                diff2 = np.nansum(np.abs(q_cur - q_old) / np.abs(q_cur)) / n_cell_sub

        p_by_t[:, gg_idx] = p_cur

        # Store q values for cells of this type
        cell_indices = np.where(sel)[0]
        q_vec[cell_indices] = q_cur

        if verbose:
            print(f"{gg} completed")
            print(f"diff1 = {diff1}")
            print(f"diff2 = {diff2}")

    return {"p_by_t": p_by_t, "q_vec": q_vec}
