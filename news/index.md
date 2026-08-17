# Changelog

## PICsnATAC 0.3.2

- Repair the chromosome-streaming `PIC_counting(load_full = FALSE)`
  path, including empty chromosomes, chromosome-local counting, and
  restoration of the input peak order.
- Make the DAR artifact threshold explicit and configurable. Counts
  strictly above `artifact_threshold` are treated as unreliable; the
  remaining counts are capped at the model-supported value of 5.
- Accept both base and sparse matrices consistently in insertion-rate
  helpers.
- Binarize dense and sparse inputs correctly in
  [`get_r_by_ct_mat_pq()`](https://zhen-miao.github.io/PICsnATAC/reference/get_r_by_ct_mat_pq.md)
  and add stable convergence and input validation.
- Expand automated tests, package metadata, and continuous integration
  checks.
- Declare the `R.utils` runtime dependency required by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  when reading gzip-compressed fragment files.
