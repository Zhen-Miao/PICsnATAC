# Get region (peak) by cell type matrix and per cell capturing rate

For a snATAC-seq (binary) dataset, compute the peak-specific open
probability and cell-specific capturing rates

## Usage

``` r
get_r_by_ct_mat_pq(
  cell_type_set,
  r_by_c,
  cell_type_labels,
  n_features_per_cell,
  p_acc = 5e-04,
  q_acc = 5e-04,
  n_max_iter = 800,
  verbose = TRUE
)
```

## Arguments

- cell_type_set:

  A vector containing all cell types

- r_by_c:

  Input matrix, region (peak) by cell

- cell_type_labels:

  A vector containing cell type labels

- n_features_per_cell:

  The number of features in the matrix, can be calculated by
  nrow(r_by_c)

- p_acc:

  The accuracy of p, default specified as 0.0005

- q_acc:

  The accuracy of q, default specified as 0.0005

- n_max_iter:

  The maximum iteration, default = 800

- verbose:

  Whether to output information on processing status

## Value

A list with two elements,

- p_by_t Peak by cell type matrix, each element represents the open
  probability of the peak in the corresponding cell type

- q_vec A vector of cell-specific capturing rate
