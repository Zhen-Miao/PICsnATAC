
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

low = readRDS('GSE162690_CellLine_LowLoading-PeakMatrix-SE.rds')
cells <- colnames(low@assays$data@listData$PeakMatrix)
cells <- gsub(".*#","",cells)
peak_sets = low@rowRanges
rm(low)

mtxsig <- CreateFragmentObject(
  path = "/test_github/GSE162690_CellLine_LowLoading.fragments.tsv.gz",
  cells = cells
)


mtxsig2 <- FeatureMatrix(
  fragments = mtxsig,
  features = peak_sets,
  cells = cells
)

saveRDS(mtxsig2, 'mtx_signac_low_archR.rds')



