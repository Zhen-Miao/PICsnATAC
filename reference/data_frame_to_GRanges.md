# convert peak_set into GRanges object if input is a data.frame

convert peak_set into GRanges object if input is a data.frame

## Usage

``` r
data_frame_to_GRanges(peak_sets)
```

## Arguments

- peak_sets:

  A data.frame object of peaks that we want to use as features The first
  column should be seqname, (e.g., chr1); the second column should be
  the start site, and the third column should be the end site. This
  should be after 5 bp and 4 bp correction of Tn5 insertion location.

## Value

A GRanges object that contain the same information as peak_sets
