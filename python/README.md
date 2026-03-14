# PICsnATAC (Python)

Python implementation of PICsnATAC: **Paired Insertion Counting for single-nucleus ATAC-seq data**.

This is a Python port of the [R package](https://github.com/Zhen-Miao/PICsnATAC) published in *Nature Methods* 21.1 (2024): 32-36.

## Installation

```bash
pip install .
```

### Dependencies

- numpy, pandas, scipy
- pysam (replaces R's Rsamtools for tabix-indexed file access)
- tqdm (progress bars)
- joblib (parallel processing)

## Quick Start

```python
import pandas as pd
from picsnATAC import PIC_counting

# Load cell barcodes
meta = pd.read_csv("atac_pbmc_5k_nextgem_singlecell.csv")
cells = meta.loc[meta["is__cell_barcode"] == 1, "barcode"].tolist()

# Load peaks
peaks = pd.read_csv("atac_pbmc_5k_nextgem_peaks.bed", sep="\t",
                     header=None, names=["chrom", "start", "end"], comment="#")

# Run PIC counting
result = PIC_counting(
    cells=cells,
    fragment_tsv_gz_file_location="atac_pbmc_5k_nextgem_fragments.tsv.gz",
    peak_sets=peaks,
)

# result.matrix is a scipy.sparse.csc_matrix (peaks x cells)
# result.peak_names is a list of "chr:start-end" strings
# result.cell_names is a list of cell barcode strings
```

## API Reference

### PIC Counting

- **`PIC_counting(cells, fragment_tsv_gz_file_location, peak_sets, ...)`** - Build peak-by-cell PIC count matrix
- **`load_fragments(fragment_tsv_gz_file_location, cells)`** - Load and filter fragment file
- **`count_peaks(...)`** - Count paired insertions for a single cell

### Parameter Estimation

- **`get_r_by_ct_mat_pq(cell_type_set, r_by_c, cell_type_labels, ...)`** - EM algorithm for peak open probability and capturing rate
- **`obs_to_insertion_MLE_lam(pic_mat, capturing_rates, plen, n_cores)`** - MLE for insertion rates
- **`obs_to_insertion_MLE_obj(pic_mat, capturing_rates, plen, n_cores)`** - Max log-likelihood values
- **`obs_to_insertion_ME(pic_mat, capturing_rates, cell_type_labels, ...)`** - Moment estimator for insertion rates

### Differential Accessibility

- **`DAR_by_LRT(pic_mat, capturing_rates, cell_type_labels, ...)`** - Likelihood ratio test for differential accessible regions

## Key Differences from R Version

| Feature | R | Python |
|---------|---|--------|
| Genomic ranges | GenomicRanges/Rsamtools | numpy binary search + pysam |
| Sparse matrices | Matrix::dgCMatrix | scipy.sparse.csc_matrix |
| Return type | Matrix with rownames/colnames | `PICResult` dataclass |
| Parallel processing | parallel::mclapply | joblib.Parallel |
| Fragment loading (tabix) | Rsamtools::scanTabix | pysam.TabixFile |

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```
