# Count number of paired insertions in each peak

Count number of paired insertions in each peak

## Usage

``` r
count_peaks(peak_sets, filtered_fragments, extend_size, n_features)
```

## Arguments

- peak_sets:

  GRanges object of peak sets

- filtered_fragments:

  filtered fragment also as a GRanges object

- extend_size:

  How long should we extend the exact insertion site as accessible
  window

- n_features:

  Number of features (peaks)

## Value

A sparse vector of PIC for each peak
