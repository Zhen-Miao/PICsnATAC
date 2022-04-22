
# library(SnapATAC)
library(viridisLite);
library(ggplot2);
library(rtracklayer)
library(dplyr)
library('GenomicRanges')
library('Matrix')
library('Signac')

## index
# tabix -p bed GSE162690_CellLine_LowLoading.fragments.tsv.gz

## first we do 28
barcodes <- readRDS('cells_sel.rds')
kidsample <- readRDS('kidsample.rds')
cells = barcodes[kidsample == '90028']


peak_sets = readRDS('peak_list_P0.rds')

mtxsig <- CreateFragmentObject(
  path = "/P0kidney/fragments_28.tsv.gz",
  cells = cells
)


mtxsig2 <- FeatureMatrix(
  fragments = mtxsig,
  features = peak_sets,
  cells = cells
)

saveRDS(mtxsig2, 'mtx_signac_P0_28.rds')

4750478
44837096

