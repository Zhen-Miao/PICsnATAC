# Calculate maximized log-likelihoods from observed PIC counts

Calculate maximized log-likelihoods from observed PIC counts

## Usage

``` r
obs_to_insertion_MLE_obj(pic_mat, capturing_rates, plen, n_cores = 1L)
```

## Arguments

- pic_mat:

  The observed peak by cell PIC count matrix

- capturing_rates:

  A vector of estimated capturing rates for each cell

- plen:

  A vector of peak widths

- n_cores:

  A positive integer specifying the number of cores. On Windows, values
  greater than one fall back to serial evaluation.

## Value

A numeric vector containing the maximized log-likelihood for each peak
over insertion rates from 0.01 to 20.
