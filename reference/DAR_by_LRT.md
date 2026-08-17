# Compute p value for DAR test between two cell types

Compute p value for DAR test between two cell types

## Usage

``` r
DAR_by_LRT(
  pic_mat,
  capturing_rates,
  cell_type_labels,
  n_cores = 1L,
  plen = NULL,
  min_frag_length = 25,
  max_frag_length = 600,
  estimation_approach = "MLE",
  artifact_threshold = 10L
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

  A positive integer specifying the number of cores. On Windows, values
  greater than one fall back to serial evaluation. Default = 1.

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

- artifact_threshold:

  A heuristic upper threshold for unreliable counts. Counts strictly
  greater than this value may reflect mapping errors or unusual fragment
  structures and are replaced by zero before the remaining counts are
  capped at 5 for the likelihood model. The default is 10; use `Inf` to
  disable artifact replacement. Because this cutoff is assay- and
  pipeline-dependent, sensitivity analyses with alternative values are
  recommended.

## Value

A numeric vector containing one likelihood-ratio-test p-value per peak,
named from `rownames(pic_mat)` when available.
