
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PICsnATAC

<!-- badges: start -->
<!-- badges: end -->

The goal of PICsnATAC is to construct cell by peak matrix with Paired
Insertion Counting (PIC) for snATAC-seq data

## Dependencies

Please install the dependent libraries by running the following codes
``` r
install.packages('data.table') ## (please make sure it is newer than 1.8)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```

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

[Vignette 2: PIC counting in Seurat/Signac workflow.](https://htmlpreview.github.io/?https://github.com/Zhen-Miao/PIC-snATAC/blob/main/vignettes/vignette-2----PIC-counting-in-Seurat-workflow.html)

[Vignette 3: PIC counting in ArchR workflow.](https://htmlpreview.github.io/?https://github.com/Zhen-Miao/PIC-snATAC/blob/main/vignettes/vignette-3----PIC-counting-in-ArchR-workflow.html)

[vignette 4: PIC counting with dsc-ATAC-seq data.](https://htmlpreview.github.io/?https://github.com/Zhen-Miao/PIC-snATAC/blob/main/vignettes/vignette-4----PIC-counting-with-dsc-ATAC-seq-data.html)

[vignette 5: DAR test with PIC parametric test framework.](https://htmlpreview.github.io/?https://github.com/Zhen-Miao/PIC-snATAC/blob/main/vignettes/vignette-5----DAR-analysis-with-PIC-parametric-framework.html)
