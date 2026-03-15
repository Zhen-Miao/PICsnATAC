# Compute p value for DAR test between two cell types

Compute p value for DAR test between two cell types

## Usage

``` r
DAR_by_LRT(
  pic_mat,
  capturing_rates,
  cell_type_labels,
  n_cores = 1,
  plen = NULL,
  min_frag_length = 25,
  max_frag_length = 600,
  estimation_approach = "MLE"
)
```

## Arguments

- pic_mat:

  The observed peak by cell PIC count matrix

- capturing_rates:

  A vector of estimated capturing rates for each cell

- cell_type_labels:

  A vector specifying cell type labels

- n_cores:

  A numerical value to specify the number of cores in parallel, default
  = 1

- plen:

  A vector of peak length

- min_frag_length:

  The value for the s1 hyperparameter in the ssPoisson distribution,
  this stands for the minimum fragment length requirement such that the
  fragment can be amplifiable and mappable to genome. Default = 25

- max_frag_length:

  The value for the s2 hyperparameter in the ssPoisson distribution,
  this stands for the max fragment length requirement such that the
  fragment can be amplifiable. Default = 600

- estimation_approach:

  The approach for parameter estimation, either 'MLE' for condition 1 or
  'ME' for condition 1+2. The 'MLE' approach is more accurate and
  usually it has a higher power, but it ignores the size filtering step
  in snATAC-seq data generation. Default is 'MLE'

## Value

Log loss value
