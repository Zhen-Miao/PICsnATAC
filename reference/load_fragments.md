# Function for loading fragment files and filter by cell barcodes

Function for loading fragment files and filter by cell barcodes

## Usage

``` r
load_fragments(fragment_tsv_gz_file_location, cells, verbose = TRUE)
```

## Arguments

- fragment_tsv_gz_file_location:

  The 10X Cell Ranger output fragment.tsv.gz file location. This can
  usually be found at the /out directory from Cell Ranger output

- cells:

  Cell barcode labels as a character vector.

- verbose:

  Whether to output progress message. Default TRUE

## Value

A data frame containing fragments filtered by cell barcode.
