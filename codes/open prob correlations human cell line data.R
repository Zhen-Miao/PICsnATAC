
## load data
library('Rfast')
# library("SnapATAC");
library(viridisLite);
library(ggplot2);
library("RColorBrewer")
library(dplyr)
library(reshape2)
library(aricode)
library('GenomicRanges')
library(reshape2)
# 

low = readRDS('GSE162690_CellLine_LowLoading-PeakMatrix-SE.rds')
load('archR_data_p_q_saved.RData')
features_to_keep <- readRDS('features_to_keep_low.rds')

## first check, for a particular cell type (PT), the number of two and number of cells accessible
#

## for all cells 
ctypes = unique(low@colData$DemuxletClassify2)
ctypes = setdiff(ctypes, c('DBL', 'AMB'))

pmatpt <- low@assays$data@listData$PeakMatrix[features_to_keep,]

pks = low@rowRanges[features_to_keep]
length(pks)


pmatptbin <- pmatpt
pmatptbin@x <- rep(1, times = length(pmatptbin@x))

## 12 vs 34
pmatpt@x[pmatpt@x <= 2] = 1
pmatpt@x[pmatpt@x >= 3] = 2
# -- end 12 vs 34

# ## 1 vs 2, note, we need to modify so that the quantity is correct
# 
# pmatpt@x[pmatpt@x >= 3] = 0
# pmatpt = drop0(pmatpt)
# pmatptbin <- pmatpt
# pmatptbin@x <- rep(1, times = length(pmatptbin@x))
# ## -- end 1 vs 2

proportion_counts_more_than_1 <- matrix(nrow = dim(pmatpt)[1], ncol = length(ctypes))
colnames(proportion_counts_more_than_1) <- ctypes

n_cell_acc_prop <- proportion_counts_more_than_1

## only look at 2 vs 1 
for(jj in ctypes){
  pmatpt1 = pmatpt[,low@colData$DemuxletClassify2 == jj]
  pmatptbin1 = pmatptbin[,low@colData$DemuxletClassify2 == jj]
  
  n_cell_acc <- Matrix::rowSums(pmatptbin1)
  n_counts_more_than_1 <- Matrix::rowSums(pmatpt1) - n_cell_acc
  proportion_counts_more_than_1[,jj] <- n_counts_more_than_1 / n_cell_acc
  n_cell_acc_prop[,jj] <- p_by_t_new[,jj]
}

corr_across_features <- vector(length = dim(pmatpt)[1])
for(i in 1:dim(pmatpt)[1]){
  corr_across_features[i] = cor(n_cell_acc_prop[i,], proportion_counts_more_than_1[i,],
                                method = 'spearman')
}

corr_p_across_features <- vector(length = dim(pmatpt)[1])
for(i in 1:dim(pmatpt)[1]){
  if(anyNA(proportion_counts_more_than_1[i,])){
    corr_p_across_features[i] = NA
    next
  }

  corr_p_across_features[i] = cor.test(n_cell_acc_prop[i,], proportion_counts_more_than_1[i,],
                                method = 'spearman')$p.value
}

# sum(is.na(corr_p_across_features))
# [1] 119643
p_no_na = corr_p_across_features[!is.na(corr_p_across_features)] ## 46499 peaks
padj = p.adjust(p_no_na, method = 'fdr')
sum(padj < 0.05) # 4366
sum(p_no_na < 0.05) # 16039

## for 1 and 2
p_no_na = corr_p_across_features[!is.na(corr_p_across_features)] ## 47644 peaks
padj = p.adjust(p_no_na, method = 'fdr')
sum(padj < 0.05) # 38 , 0.0007975821
sum(p_no_na < 0.05) # 5096


## for 1 vs 2 
quantile(corr_across_features,na.rm = T)
# 0%         25%         50%         75%        100% 
# -0.98480698 -0.34042710 -0.06154575  0.22424242  0.99696509 
sum(corr_across_features < 0 , na.rm = T)  # 33005
sum(corr_across_features > 0 , na.rm = T) ## 26296
quantile(corr_p_across_features,na.rm = T)

## for 3 vs 4
quantile(corr_across_features,na.rm = T)
# 0%        25%        50%        75%       100% 
# -0.7294867  0.3742613  0.5276493  0.7006490  0.9969651
sum(corr_across_features < 0 , na.rm = T)  # 2535
sum(corr_across_features > 0 , na.rm = T) ## 53169

## colors 

#f2a285 rose 
#8491b4 dark purple
#00a087 green
#3c5488 dark blue
#4dbbd5 light blue 
#00a087 green


hist(corr_across_features,xlim = c(-1,1))
hist(corr_p_across_features, xlim = c(0,1))


df = data.frame(corr = corr_across_features)

p = ggplot(df, aes(x = corr)) + 
  geom_histogram(bins = 20, color='black', fill="#8491b4") + 
  xlim(-1,1)+
  # scale_x_discrete(drop = FALSE, limits = c(-1,1))+
  ggtitle("Correlation ( p, freq(y=2|y=1 or 2) )") +
  xlab("Spearman Correlation") +
  ylab('Frequency')+
  theme_Pub()+
  scale_colour_Publication()
plot(p)

pdf('frequency of correlation f2 over f12 2.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()


p_plot = corr_p_across_features[!is.na(corr_p_across_features)]


pdf('qq_plot_f2_over_f12.pdf', height = 4, width = 4)
qqunif.plot(p_plot[p_plot != 0])
dev.off()

p = ggplot(df, aes(x = corr)) + 
  geom_histogram(bins = 20, color='black', fill="#4dbbd5") + 
  xlim(-1,1)+
  # scale_x_discrete(drop = FALSE, limits = c(-1,1))+
  ggtitle("Correlation ( p, freq(y>=3|y>0) )") +
  xlab("Spearman Correlation") +
  ylab('Frequency')+
  theme_Pub()
plot(p)

pdf('frequency of correlation 12 vs 34_new 2.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()

## example of 


df2 <- data.frame(open_prob = n_cell_acc_prop[46157,], proportion_more_than_1 = proportion_counts_more_than_1[46157,],
                  cell_type = as.factor(colnames(n_cell_acc_prop)))
p = ggplot(df2) + 
  geom_point(aes(x = open_prob, y = proportion_more_than_1, colour = cell_type), 
             size = 2, alpha = 1) + 
  xlim(0.25,0.75)+
  ylim(0.5,1)+
  # scale_x_discrete(drop = FALSE, limits = c(-1,1))+
  ggtitle("chr14:54863383-54863883") +
  xlab("Open probability in each cell type") +
  ylab('freq(y=2|y=1 or 2)')+
  theme_Pub()
plot(p)

pdf('example of frequency of 1 and open prob_46157.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()


df2 <- data.frame(open_prob = n_cell_acc_prop[107250,], proportion_more_than_1 = proportion_counts_more_than_1[107250,],
                  cell_type = as.factor(colnames(n_cell_acc_prop)))
p = ggplot(df2) + 
  geom_point(aes(x = open_prob, y = proportion_more_than_1, colour = cell_type), 
             size = 2, alpha = 1) + 
  xlim(0.1,0.8)+
  ylim(0,0.25)+
  # scale_x_discrete(drop = FALSE, limits = c(-1,1))+
  ggtitle("chr3:156272734-156273234") +
  xlab("Open probability in each cell type") +
  ylab('freq(y>=3|y>0)')+
  theme_Pub()
plot(p)

pdf('example of frequency of 34 and open prob_107250.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()


pks[107250]





