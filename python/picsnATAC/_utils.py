"""Genomic interval helper functions for PICsnATAC.

Replaces R's GenomicRanges/IRanges with numpy-based binary search operations.
"""

import numpy as np
import pandas as pd


def validate_peaks(peaks_df):
    """Validate and normalize a peaks DataFrame.

    Ensures columns are named 'chrom', 'start', 'end'.
    Accepts either those names or the first 3 columns positionally.

    Parameters
    ----------
    peaks_df : pd.DataFrame
        DataFrame with genomic peak coordinates.

    Returns
    -------
    pd.DataFrame
        Copy with standardized column names.
    """
    peaks = peaks_df.copy()
    required = {"chrom", "start", "end"}
    # Also accept 'seqname' / 'seqnames' as chrom column
    col_map = {}
    for c in peaks.columns:
        cl = c.lower()
        if cl in ("seqname", "seqnames", "chr", "chrom", "chromosome"):
            col_map[c] = "chrom"
        elif cl == "start":
            col_map[c] = "start"
        elif cl == "end":
            col_map[c] = "end"

    if required.issubset(set(col_map.values())):
        peaks = peaks.rename(columns=col_map)
    else:
        # Fall back to positional
        if peaks.shape[1] < 3:
            raise ValueError("peaks_df must have at least 3 columns")
        cols = list(peaks.columns)
        peaks = peaks.rename(columns={cols[0]: "chrom", cols[1]: "start", cols[2]: "end"})

    # Drop rows with missing values (e.g., from comment lines)
    peaks = peaks.dropna(subset=["chrom", "start", "end"])
    peaks["start"] = peaks["start"].astype(np.int64)
    peaks["end"] = peaks["end"].astype(np.int64)
    peaks["chrom"] = peaks["chrom"].astype(str)

    return peaks[["chrom", "start", "end"]].reset_index(drop=True)


def find_overlaps_first(query_chrom, query_start, query_end,
                        subject_chrom, subject_start, subject_end):
    """For each query interval, find the index of the first overlapping subject.

    Uses binary search (np.searchsorted) per chromosome.
    Peaks (subjects) must be sorted within each chromosome.

    Parameters
    ----------
    query_chrom, query_start, query_end : np.ndarray
        Query intervals.
    subject_chrom, subject_start, subject_end : np.ndarray
        Subject (peak) intervals, sorted by start within each chromosome.

    Returns
    -------
    np.ndarray
        Integer array of length len(query). Contains the global subject index
        for the first overlap, or -1 if no overlap.
    """
    n_query = len(query_chrom)
    result = np.full(n_query, -1, dtype=np.int64)

    # Group subjects by chromosome
    unique_chroms = np.unique(subject_chrom)
    subj_chrom_indices = {}
    for chrom in unique_chroms:
        mask = subject_chrom == chrom
        indices = np.where(mask)[0]
        subj_chrom_indices[chrom] = indices

    # Process queries per chromosome
    for chrom in unique_chroms:
        if chrom not in subj_chrom_indices:
            continue
        subj_idx = subj_chrom_indices[chrom]
        s_starts = subject_start[subj_idx]
        s_ends = subject_end[subj_idx]

        q_mask = query_chrom == chrom
        if not np.any(q_mask):
            continue
        q_indices = np.where(q_mask)[0]
        q_starts = query_start[q_indices]
        q_ends = query_end[q_indices]

        # For overlap: query_start < subject_end AND query_end > subject_start
        # Find first subject whose end > query_start
        pos = np.searchsorted(s_ends, q_starts, side="right")

        # Check bounds and verify overlap
        valid = pos < len(subj_idx)
        valid_idx = np.where(valid)[0]
        if len(valid_idx) > 0:
            candidate_starts = s_starts[pos[valid_idx]]
            overlaps = q_ends[valid_idx] > candidate_starts
            hit_idx = valid_idx[overlaps]
            result[q_indices[hit_idx]] = subj_idx[pos[hit_idx]]

    return result


def subset_by_overlaps_mask(frag_chrom, frag_start, frag_end,
                            peak_chrom, peak_start, peak_end,
                            maxgap=0):
    """Return boolean mask for fragments overlapping any peak.

    Parameters
    ----------
    frag_chrom, frag_start, frag_end : np.ndarray
        Fragment coordinates.
    peak_chrom, peak_start, peak_end : np.ndarray
        Peak coordinates (sorted by start within each chromosome).
    maxgap : int
        Maximum gap between fragment and peak to still count as overlap.

    Returns
    -------
    np.ndarray
        Boolean mask of length len(frag_chrom).
    """
    # Expand peaks by maxgap for the overlap test
    expanded_start = peak_start - maxgap
    expanded_end = peak_end + maxgap
    hits = find_overlaps_first(
        frag_chrom, frag_start, frag_end,
        peak_chrom, expanded_start, expanded_end,
    )
    return hits >= 0


def parse_peak_names(peak_names):
    """Parse 'chr1:1000-2000' formatted peak names.

    Parameters
    ----------
    peak_names : array-like of str
        Peak name strings.

    Returns
    -------
    tuple of (chrom, start, end) as numpy arrays
    """
    chroms = []
    starts = []
    ends = []
    for name in peak_names:
        chrom, rest = name.split(":")
        start, end = rest.split("-")
        chroms.append(chrom)
        starts.append(int(start))
        ends.append(int(end))
    return np.array(chroms), np.array(starts, dtype=np.int64), np.array(ends, dtype=np.int64)
