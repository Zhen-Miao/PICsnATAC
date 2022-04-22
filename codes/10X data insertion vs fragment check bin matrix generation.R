## process 10X data


require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
# library(Signac)
library('ggplot2')
library(ggsci)

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

## analysis with ArchR
library('ArchR')
library(parallel)
addArchRThreads(threads = 4) 
addArchRGenome("hg38")

inputFiles <- 'atac_pbmc_5k_nextgem_fragments.tsv.gz'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = 'pbmc_5k',
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 500, 
  addTileMat = F,
  addGeneScoreMat = F
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/pbmc_atac_arrow2",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

projHeme2 <- addTileMatrix(input = projHeme1,
                           binarize = F, force = T)

# str(projHeme2)
# 
# 
# projHeme2@cellColData
# 
# getSampleColData(projHeme2)[,"ArrowFiles"]

ss = getMatrixFromProject(projHeme2,useMatrix = 'TileMatrix')

str(ss)
bins = ss@elementMetadata@listData
bindf = data.frame(sname = base::rep(bins$seqnames@values, times = bins$seqnames@lengths),
                   start = bins$start)
tm = ss@assays@data$TileMatrix

cells.sel = paste('pbmc_5k#', colnames(mtx), sep = ''  )


cs = intersect(cells.sel, colnames(tm))
tm = tm[,cs]

rs = rowSums(tm)

quantile(rs)
tm = tm[rs>=10,]
bindf = bindf[rs >= 10,]

table(tm@x)

saveRDS(tm, 'tile matrix filtered.rds')

ta <- table(tm@x[tm@x <= 8])
df1 = as.data.frame(ta)
labs = format(round(df1$Freq /(10^7), 2),nsmall = 2)

p = ggplot(df1[1:8,], aes(x=Var1, y = Freq/(10^7))) +
  geom_bar(stat="identity", color="darkblue", fill="#00A087B2") + 
  geom_text(aes(label=labs, vjust=-0.2, hjust = 0.5), size=3)+
  xlab("Count") +
  ylab('Frequency (X10^7)')+
  ggtitle("500bp bin matrix, insertion-based") +
  theme_Pub()

plot(p)


pdf('frequency count bin ArchR.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()


cells = gsub(".*#","",cs)
saveRDS(cells,'cells_selected_bin_mat.rds')
cells = readRDS('cells_selected_bin_mat.rds')

library('Signac')
mtxsig <- CreateFragmentObject(
  path = "/10X_PBMC_ATAC/atac_pbmc_5k_nextgem_fragments.tsv.gz",
  cells = cells
)

# bindf$end = bindf$start + 499
# colnames(bindf) = c('seqname','start','end')
# bin_sets = GenomicRanges::makeGRangesFromDataFrame(bindf)
# 
# saveRDS(bin_sets, 'bin_sets_pbmc5k_ATAC.rds')
bin_sets = readRDS('bin_sets_pbmc5k_ATAC.rds')

mtxsig3 <- FeatureMatrix(
  fragments = mtxsig,
  features = bin_sets,
  cells = cells
)

saveRDS(mtxsig3, 'mtx_from_signac_bin.rds')
mtxsig3 = readRDS('mtx_from_signac_bin.rds')

ta2 <- table(mtxsig3@x[mtxsig3@x <= 8])
df1 = as.data.frame(ta2)
labs = format(round(df1$Freq /(10^7), 2),nsmall = 2)
labs[6:8] <- '<0.01'

p = ggplot(df1[1:8,], aes(x=Var1, y = Freq/(10^7))) +
  geom_bar(stat="identity", color="darkblue", fill="#00A087B2") + 
  geom_text(aes(label=labs, vjust=-0.2, hjust = 0.5), size=3)+
  xlab("Count") +
  ylab('Frequency (X10^7)')+
  ggtitle("500bp bin matrix, fragment-based") +
  theme_Pub()

plot(p)


pdf('frequency count bin signac 2.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()





# 
# projHeme1
# getAvailableMatrices(projHeme1)
# 
# projHeme1[projHeme1$cellNames[1:100], ]
# 
# str(projHeme1@cellColData)
# projHeme1
# ArchR::addTileMatrix()


# mtxsig <- CreateFragmentObject(
#   path = "/Users/zhenmiao/Dropbox/subspace clustering visualization project/10X_PBMC_ATAC/atac_pbmc_5k_nextgem_fragments.tsv.gz",
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


mypal = pal_aaas("default", alpha = 0.7)(9)
mypal
library("scales")
show_col(mypal)

ta <- table(mtx@x[mtx@x <= 8])
df1 = as.data.frame(ta)

labs = format(round(df1$Freq /(10^7), 2),nsmall = 2)

p = ggplot(df1[1:8,], aes(x=Var1, y = Freq/(10^7))) +
  geom_bar(stat="identity", color="#5F559BB2", fill="#008B45B2") + 
  geom_text(aes(label=labs, vjust=-0.2, hjust = 0.5), size=3)+
  xlab("Count") +
  ylab('Frequency (X10^7)')+
  ggtitle("peak matrix, insertion-based") +
  theme_Publication()+
  scale_colour_aaas()

plot(p)


pdf('frequency count peak insertion 2 .pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()



ta2 <- table(mtxsig2@x[mtxsig2@x <= 8])
df2 = as.data.frame(ta2)

labs = format(round(df2$Freq /(10^7), 2),nsmall = 2)
labs[7:8] <- '<0.01'

p = ggplot(df2[1:8,], aes(x=Var1, y = Freq/(10^7))) +
  geom_bar(stat="identity", color="#5F559BB2", fill="#008B45B2") + 
  geom_text(aes(label=labs, vjust=-0.2, hjust = 0.5), size=3)+
  xlab("Count") +
  ylab('Frequency (X10^7)')+
  ggtitle("peak matrix, fragment-based") +
  theme_Publication()+
  scale_colour_aaas()

plot(p)


pdf('frequency count peak fragment 2.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()




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
  theme(panel.background = element_blank(),
        plot.background = element_rect(colour = NA),
        # panel.border = element_rect(colour = NA),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        axis.title = element_text(face = "bold",size = rel(1)),
        legend.key = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.size= unit(0.4, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(face="italic"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"))+
  scale_colour_Publication()

plot(p)

pdf('frequency count of Signac.pdf', width = 5, height = 5)
plot(p)
dev.off()


## bquote('Frequency' ~(X10^7))










