## process 10X data


require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
library(Signac)
library('ggplot2')

mex_dir_path <- "/10X_PBMC_ATAC/filtered_peak_bc_matrix"


mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

dim(mtx)


# mtxsig <- CreateFragmentObject(
#   path = "/10X_PBMC_ATAC/atac_pbmc_5k_nextgem_fragments.tsv.gz",
#   cells = barcodes$barcode
# )
# 
# feature_df = read.table('./filtered_peak_bc_matrix/peaks.bed')
# colnames(feature_df) <- c('seqname', 'start', 'end')
# peak_sets = GenomicRanges::makeGRangesFromDataFrame(feature_df)
# median(GenomicRanges::width(peak_sets))
# 
# mtxsig2 <- FeatureMatrix(
#   fragments = mtxsig,
#   features = peak_sets,
#   cells = barcodes$barcode
# )
# 
# saveRDS(mtxsig2, 'mtx_from_signac.rds')

mtxsig2 <- readRDS('mtx_from_signac.rds')


ta <- table(mtx@x[mtx@x <= 10])
df1 = as.data.frame(ta)



ta2 <- table(mtxsig2@x[mtxsig2@x <= 10])
df2 = as.data.frame(ta2)

plot(as.numeric(df2$Var1), log(df2$Freq))


df21 = df2[1:6,]
df21$Freq = df21$Freq / (10^7)

log.model <-lm(log(Freq) ~ as.numeric(Var1), df21)
log.model.df <- data.frame(x = as.numeric(df21$Var1),
                           y = exp(fitted(log.model)))


## check some quantities 
p = ggplot(df2[1:6,], aes(x=Var1, y = Freq/(10^7))) +
  geom_bar(stat="identity", color="darkblue", fill="lightblue") + 
  geom_line(data = log.model.df, aes(x, y, color = "Log Model"),
            size = 1, linetype = 1) + 
  geom_text(aes(label=Freq), vjust=-0.3, size=3)+
  xlab("Count") +
  ylab('Frequency (X10^7)')+
  theme_Pub()

plot(p)

pdf('frequency count of Signac.pdf', width = 5, height = 5)
plot(p)
dev.off()


## bquote('Frequency' ~(X10^7))










