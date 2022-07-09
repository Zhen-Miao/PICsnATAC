
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PICsnATAC

<!-- badges: start -->
<!-- badges: end -->

The goal of PICsnATAC is to construct cell by peak matrix with Paired
Insertion Counting (PIC) for snATAC-seq data

## Installation

You can install the development version of PICsnATAC from
[GitHub](https://github.com/Zhen-Miao/PIC-snATAC) with:

``` r
# install.packages("devtools")
devtools::install_github("Zhen-Miao/PIC-snATAC")
```

## Example

This is a basic example which shows you how to construct PIC-based cell
by peak matrix:

``` r
library('PICsnATAC')

pic_matrix <- PIC_counting(cells, fragment_tsv_gz_file_location, peak_sets)
```
For additional example, please check the following vignettes 
[Vignette 1: PIC counting for Cell Ranger output.](https://htmlpreview.github.io/?https://github.com/Zhen-Miao/PIC-snATAC/blob/main/vignettes/vignette-1----PIC-counting-with-10X-Cell-Ranger-output.html)
