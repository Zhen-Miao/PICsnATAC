# Moment estimator for insertion rates from observed values

Moment estimator for insertion rates from observed values

## Usage

``` r
obs_to_insertion_ME(
  pic_mat,
  capturing_rates,
  cell_type_labels,
  cap_insertion = 20,
  min_frag_length = 25,
  max_frag_length = 600,
  insertion_rates = (1:2000) * 0.01,
  peak_lengths = 4:20 * 50
)
```

## Arguments

- pic_mat:

  The observed peak by cell PIC count matrix

- capturing_rates:

  A vector of estimated capturing rates for each cell

- cell_type_labels:

  A vector of cell type lables for each cell

- cap_insertion:

  The maximum number of insertions in a peak region.

- min_frag_length:

  The value for the s1 hyperparameter in the ssPoisson distribution,
  this stands for the minimum fragment length requirement such that the
  fragment can be amplifiable and mappable to genome. Default = 25

- max_frag_length:

  The value for the s2 hyperparameter in the ssPoisson distribution,
  this stands for the max fragment length requirement such that the
  fragment can be amplifiable. Default = 600

- insertion_rates:

  The range of insertion rates (per 1000 bp) to be considered.

- peak_lengths:

  A vector of peak lengths

## Value

A matrix of estimated insertion rates
