"""Core PIC counting functions for snATAC-seq data.

Translates R/PIC_function.R to Python.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pysam
from scipy import sparse
from tqdm import tqdm

from picsnATAC._utils import (
    find_overlaps_first,
    subset_by_overlaps_mask,
    validate_peaks,
)


@dataclass
class PICResult:
    """Result container for PIC_counting, since scipy sparse matrices lack row/col names."""
    matrix: sparse.csc_matrix
    peak_names: list
    cell_names: list


def load_fragments(fragment_tsv_gz_file_location, cells, verbose=True):
    """Load fragment file and filter by cell barcodes.

    Parameters
    ----------
    fragment_tsv_gz_file_location : str
        Path to the fragment.tsv.gz file.
    cells : list of str
        Cell barcode strings to retain.
    verbose : bool
        Whether to print progress messages.

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, cell_barcode
    """
    f1 = pd.read_csv(
        fragment_tsv_gz_file_location,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chrom", "start", "end", "cell_barcode"],
        comment="#",
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "cell_barcode": str},
    )

    cell_set = set(cells)
    mask = f1["cell_barcode"].isin(cell_set)
    n_cells_found = mask.sum()

    if verbose:
        prop = n_cells_found / len(f1)
        print(f"proportion of reads in cell barcodes is {prop:.2f}")

    if n_cells_found < 1:
        raise ValueError("Cell barcodes not found in fragment files, please check input")
    elif n_cells_found < 10:
        warnings.warn("Fewer than 10 cells found in fragment files, please consider checking input")

    return f1.loc[mask].reset_index(drop=True)


def count_peaks(peak_chrom, peak_start, peak_end,
                frag_chrom, frag_start, frag_end,
                extend_size, n_features):
    """Count paired insertions for one cell's fragments against peaks.

    Parameters
    ----------
    peak_chrom, peak_start, peak_end : np.ndarray
        Peak coordinates.
    frag_chrom : np.ndarray
        Fragment chromosome.
    frag_start, frag_end : np.ndarray
        Fragment start and end positions (0-based, end-exclusive BED format).
    extend_size : int
        Extension around insertion site.
    n_features : int
        Total number of peaks.

    Returns
    -------
    tuple of (indices, counts)
        Nonzero peak indices and their PIC counts.
    """
    if len(frag_chrom) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    half_ext = extend_size // 2
    half_ext_ceil = (extend_size + 1) // 2

    # Start insertion: the start position of the fragment
    # In R: resize(width=1, fix="start") gives 1bp at start, then extend centered
    # BED 0-based: insertion point is at frag_start
    ins_s_start = frag_start - half_ext
    ins_s_end = frag_start + half_ext_ceil

    # End insertion: the end position of the fragment
    # In R: resize(width=1, fix="end") gives 1bp at end position
    # R GRanges uses 1-based inclusive, so end insertion = end coord
    # In BED 0-based, the fragment end is exclusive, so insertion = frag_end - 1
    # But R reads BED 0-based start as GRanges 1-based start, and BED end as GRanges end
    # resize(width=1, fix="end") on GRanges [start, end] gives [end, end]
    # which in 0-based is [end-1, end), i.e., the position at (end-1)
    # However, since fread reads BED and then makeGRangesFromDataFrame converts:
    #   BED start -> GRanges start (adds 1 implicitly? No - GRanges just takes the value)
    # Actually: R's makeGRangesFromDataFrame just uses the values as-is for start/end
    # The fragment file has 0-based start and 0-based open end
    # R reads these as integers and creates GRanges with those exact values
    # resize(width=1, fix="end") on range [s, e] gives [e, e] (1bp at end)
    # Then resize(width=extend_size, fix="center") centers around that point
    # So end insertion point = frag_end (the value in the file)
    ins_e_start = frag_end - half_ext
    ins_e_end = frag_end + half_ext_ceil

    # Find overlapping peaks
    overlapped_s = find_overlaps_first(
        frag_chrom, ins_s_start, ins_s_end,
        peak_chrom, peak_start, peak_end,
    )
    overlapped_e = find_overlaps_first(
        frag_chrom, ins_e_start, ins_e_end,
        peak_chrom, peak_start, peak_end,
    )

    # Where both insertions overlap the same peak, set end overlap to -1
    # to avoid double-counting (this is the PIC deduplication step)
    overlapped_e_nodc = overlapped_e.copy()
    same_peak = (overlapped_s >= 0) & (overlapped_s == overlapped_e_nodc)
    overlapped_e_nodc[same_peak] = -1

    # Combine all overlaps
    ol = np.concatenate([overlapped_s, overlapped_e_nodc])
    ol = ol[ol >= 0]

    if len(ol) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    # Tabulate counts per peak
    indices, counts = np.unique(ol, return_counts=True)
    return indices, counts


def PIC_counting(cells, fragment_tsv_gz_file_location, peak_sets,
                 deduplicate=False, load_full=True, extend_size=5,
                 verbose=True):
    """Construct peak-by-cell PIC count matrix.

    Parameters
    ----------
    cells : list of str
        Cell barcode strings.
    fragment_tsv_gz_file_location : str
        Path to fragments.tsv.gz file.
    peak_sets : pd.DataFrame
        Peaks with columns chrom/start/end (or first 3 columns).
    deduplicate : bool
        Remove duplicate fragments within the same cell.
    load_full : bool
        If True, load entire fragment file into memory.
        If False, use pysam TabixFile for per-chromosome loading.
    extend_size : int
        Extension around insertion site (default 5).
    verbose : bool
        Show progress information.

    Returns
    -------
    PICResult
        Contains the sparse matrix, peak names, and cell names.
    """
    if extend_size < 0:
        raise ValueError("extend_size has to be a positive integer!")

    # Validate peaks
    peaks = validate_peaks(peak_sets)
    n_cells = len(cells)
    n_features = len(peaks)

    if n_cells == 0 or any(c is None for c in cells):
        raise ValueError("cell names are empty or contain NA values!")

    peak_chrom = peaks["chrom"].values
    peak_start = peaks["start"].values
    peak_end = peaks["end"].values

    # Generate peak name strings
    peak_names = [
        f"{c}:{s}-{e}"
        for c, s, e in zip(peak_chrom, peak_start, peak_end)
    ]

    if load_full:
        return _pic_counting_load_full(
            cells, fragment_tsv_gz_file_location, peaks,
            peak_chrom, peak_start, peak_end, peak_names,
            n_cells, n_features, deduplicate, extend_size, verbose,
        )
    else:
        return _pic_counting_tabix(
            cells, fragment_tsv_gz_file_location, peaks,
            peak_chrom, peak_start, peak_end, peak_names,
            n_cells, n_features, deduplicate, extend_size, verbose,
        )


def _pic_counting_load_full(cells, fragment_file, peaks,
                            peak_chrom, peak_start, peak_end, peak_names,
                            n_cells, n_features, deduplicate, extend_size,
                            verbose):
    """PIC counting with full file loading."""
    # Load fragments
    f1 = load_fragments(fragment_file, cells, verbose=verbose)

    # Swap start/end if start >= end (s3-ATAC-seq handling)
    bad = f1["start"] - 1 >= f1["end"]
    if bad.any():
        tmp = f1.loc[bad, "start"].values.copy()
        f1.loc[bad, "start"] = f1.loc[bad, "end"]
        f1.loc[bad, "end"] = tmp

    # Pre-filter: keep only fragments that overlap peaks
    frag_chrom = f1["chrom"].values
    frag_start = f1["start"].values
    frag_end = f1["end"].values

    overlap_mask = subset_by_overlaps_mask(
        frag_chrom, frag_start, frag_end,
        peak_chrom, peak_start, peak_end,
        maxgap=int(np.ceil(extend_size / 2)),
    )
    f1 = f1.loc[overlap_mask].reset_index(drop=True)

    if verbose:
        print("Computing peak vector for each cell.")

    # Build sparse matrix using COO format
    all_rows = []
    all_cols = []
    all_data = []

    cell_set = set(cells)
    cell_to_idx = {c: i for i, c in enumerate(cells)}

    iterator = tqdm(cells, disable=not verbose, ncols=60,
                    bar_format="[{bar}] {percentage:.0f}% finished, elapsed: {elapsed}")

    for cell in iterator:
        cell_frags = f1.loc[f1["cell_barcode"] == cell]
        if len(cell_frags) == 0:
            continue

        fc = cell_frags["chrom"].values
        fs = cell_frags["start"].values
        fe = cell_frags["end"].values

        if deduplicate:
            # Remove duplicate fragments (same chrom, start, end)
            uniq_mask = ~cell_frags.duplicated(subset=["chrom", "start", "end"])
            fc = fc[uniq_mask.values]
            fs = fs[uniq_mask.values]
            fe = fe[uniq_mask.values]

        indices, counts = count_peaks(
            peak_chrom, peak_start, peak_end,
            fc, fs, fe,
            extend_size, n_features,
        )

        if len(indices) > 0:
            col_idx = cell_to_idx[cell]
            all_rows.append(indices)
            all_cols.append(np.full(len(indices), col_idx, dtype=np.int64))
            all_data.append(counts)

    if verbose:
        print("Summarizing cell-by-peak matrix")

    # Build sparse matrix
    if len(all_rows) > 0:
        rows = np.concatenate(all_rows)
        cols = np.concatenate(all_cols)
        data = np.concatenate(all_data)
        mat = sparse.coo_matrix(
            (data, (rows, cols)),
            shape=(n_features, n_cells),
        ).tocsc()
    else:
        mat = sparse.csc_matrix((n_features, n_cells))

    return PICResult(matrix=mat, peak_names=peak_names, cell_names=list(cells))


def _pic_counting_tabix(cells, fragment_file, peaks,
                        peak_chrom, peak_start, peak_end, peak_names,
                        n_cells, n_features, deduplicate, extend_size,
                        verbose):
    """PIC counting with per-chromosome TabixFile loading."""
    tbx = pysam.TabixFile(fragment_file)

    if verbose:
        print("Data loaded by chromosome")

    cell_set = set(cells)
    cell_to_idx = {c: i for i, c in enumerate(cells)}

    # Group peaks by chromosome
    unique_chroms = []
    seen = set()
    for c in peak_chrom:
        if c not in seen:
            unique_chroms.append(c)
            seen.add(c)

    all_rows = []
    all_cols = []
    all_data = []

    chrom_iterator = tqdm(unique_chroms, disable=not verbose, ncols=60,
                          bar_format="[{bar}] {percentage:.0f}% computed, elapsed: {elapsed}")

    for chrom in chrom_iterator:
        chrom_mask = peak_chrom == chrom
        chrom_peak_start = peak_start[chrom_mask]
        chrom_peak_end = peak_end[chrom_mask]
        # Global indices for this chromosome's peaks
        chrom_peak_global_idx = np.where(chrom_mask)[0]
        chrom_peak_chrom = peak_chrom[chrom_mask]

        region_start = int(chrom_peak_start.min())
        region_end = int(chrom_peak_end.max())

        # Fetch fragment lines for this chromosome region
        try:
            lines = list(tbx.fetch(chrom, region_start, region_end))
        except ValueError:
            continue

        if not lines:
            continue

        # Parse fragment lines
        frags = []
        for line in lines:
            parts = line.split("\t")
            barcode = parts[3]
            if barcode in cell_set:
                frags.append((parts[0], int(parts[1]), int(parts[2]), barcode))

        if not frags:
            continue

        frag_df = pd.DataFrame(frags, columns=["chrom", "start", "end", "cell_barcode"])

        if verbose:
            print("Computing peak vector for each cell.")

        for cell in cells:
            cell_frags = frag_df.loc[frag_df["cell_barcode"] == cell]
            if len(cell_frags) == 0:
                continue

            fc = cell_frags["chrom"].values
            fs = cell_frags["start"].values
            fe = cell_frags["end"].values

            if deduplicate:
                uniq_mask = ~cell_frags.duplicated(subset=["chrom", "start", "end"])
                fc = fc[uniq_mask.values]
                fs = fs[uniq_mask.values]
                fe = fe[uniq_mask.values]

            # Count against ALL peaks (not just this chrom's) to match R behavior
            # where count_peaks uses the full peak_sets GRanges
            indices, counts = count_peaks(
                peak_chrom, peak_start, peak_end,
                fc, fs, fe,
                extend_size, n_features,
            )

            if len(indices) > 0:
                col_idx = cell_to_idx[cell]
                all_rows.append(indices)
                all_cols.append(np.full(len(indices), col_idx, dtype=np.int64))
                all_data.append(counts)

    if verbose:
        print("Summarizing cell-by-peak matrix")

    tbx.close()

    if len(all_rows) > 0:
        rows = np.concatenate(all_rows)
        cols = np.concatenate(all_cols)
        data = np.concatenate(all_data)
        mat = sparse.coo_matrix(
            (data, (rows, cols)),
            shape=(n_features, n_cells),
        ).tocsc()
    else:
        mat = sparse.csc_matrix((n_features, n_cells))

    return PICResult(matrix=mat, peak_names=peak_names, cell_names=list(cells))
