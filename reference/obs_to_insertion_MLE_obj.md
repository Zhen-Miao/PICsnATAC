# Calculate the optimal loss from observed PIC counts

Calculate the optimal loss from observed PIC counts

## Usage

``` r
obs_to_insertion_MLE_obj(pic_mat, capturing_rates, plen, n_cores)
```

## Arguments

- pic_mat:

  The observed peak by cell PIC count matrix

- capturing_rates:

  A vector of estimated capturing rates for each cell

- plen:

  A vector of peak widths

- n_cores:

  A numerical value to specify the number of cores in parallel

## Value

The optimized loss over insertion rates from 0.01 to 20
