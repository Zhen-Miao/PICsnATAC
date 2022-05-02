
library(viridisLite);
library(ggplot2);
library(rtracklayer)
library(dplyr)
library('GenomicRanges')
library('Matrix')
library("ggsci")
library("gridExtra")



Xc_sub = readRDS('Xc_sub_1234PBMC_sorted.rds')
Xr_sub = readRDS('Xr_sub_1234PBMC_sorted.rds')
count_sum = readRDS('count_sum_peak_gene_PBMC_sorted_1234 g3 global_filtering .rds')
peak_gene_pairs = readRDS('peak_gene_pairs_PBMC_sorted_1234 g3 global_filtering  .rds')

quantile(count_sum[,'p1_2'], na.rm = T)
quantile(count_sum[,'p12_34'], na.rm = T)

sum(count_sum[,'p1_2'] < 0.05, na.rm = T) # 463
sum(count_sum[,'p12_34'] < 0.05, na.rm = T) # 592

sum(count_sum[,'nz_prop_1'] < count_sum[,'nz_prop_2'] , na.rm = T) # 2574
sum(count_sum[,'nz_prop_1'] > count_sum[,'nz_prop_2'] , na.rm = T) # 3051

sum(count_sum[,'nz_prop_1'] > count_sum[,'nz_prop_102'] , na.rm = T) # 2877
sum(count_sum[,'nz_prop_1'] < count_sum[,'nz_prop_102'] , na.rm = T) # 2749

cdf = as.data.frame(count_sum)
cdf$fc1_2 = log(cdf$mean_2 / cdf$mean_1)
cdf$fc12_34 = log(cdf$mean_34 / cdf$mean_12)


hist(cdf[,'fc1_2'])
hist(cdf[,'fc12_34'])



countdf = cdf ## 9240   12
countdf = na.omit(countdf) ##  6931   12
countdf$fc1_2[countdf$fc1_2 == Inf | countdf$fc1_2 == -Inf ] = NA
countdf$fc12_34[countdf$fc12_34 == Inf | countdf$fc12_34 == -Inf ] = NA

## remove any row containing any NA
countdf = na.omit(countdf) ## 6619   12
quantile(countdf$nz_prop_g0)
sum(countdf$nz_prop_g0 >= 0.1) # 3387
countdf = countdf[countdf$nz_prop_g0 >= 0.1,]

## p value correction
countdf$p1_2 = p.adjust(countdf$p1_2,method = 'fdr')
countdf$p12_34 = p.adjust(countdf$p12_34,method = 'fdr')



sum(countdf$nz_prop_1 > countdf$nz_prop_2) # 1210
sum(countdf$nz_prop_1 < countdf$nz_prop_2) # 1324

sum(countdf$nz_prop_12 > countdf$nz_prop_34) # 831
sum(countdf$nz_prop_12 < countdf$nz_prop_34) # 1703


sum(countdf$p12_34 < 0.05 & countdf$fc12_34 > 0) # 189
sum(countdf$p12_34 < 0.05 & countdf$fc12_34 < 0) # 9

sum(countdf$p1_2 < 0.05 & countdf$fc1_2 > 0) # 9
sum(countdf$p1_2 < 0.05 & countdf$fc1_2 < 0) # 9



sum(countdf$p1_2 < 0.05 & countdf$fc1_2 > 0 )
sum(countdf$p1_2 < 0.05 & countdf$fc1_2 < 0 )


p1 <- ggplot(data=countdf, aes(x=fc1_2, y=-log10(p1_2))) + 
  geom_point(alpha = 1, 
             shape = 16,
             size = 1, colour = '#00A087B2') +
  annotate("text", x = 1.6, y = 14, label = "9 sig. & logFC>0", 
           family= theme_get()$text[["family"]], )+
  annotate("text", x = -1.6, y = 14, label = "9 sig. & logFC<0",
           family= theme_get()$text[["family"]], )+
  xlab("log fold change") +
  ylab('-log10(FDR corrected p value)')+
  ggtitle('Expression associated with y=1 vs y=2')+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = c(seq(-3, 3, 1)),     
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(seq(0, 15, 5)),     
                     limits = c(0, 15))   +
  theme_Pub()


pdf('volcano p vs LFC 1_2 PBMC_sort global filtering .pdf', width = 3.5, height = 3.5)
plot(p1)
dev.off()

sum(countdf$p12_34 < 0.05 & countdf$fc12_34 > 0 )
sum(countdf$p12_34 < 0.05 & countdf$fc12_34 < 0 )

p1 <- ggplot(data=countdf, aes(x=fc12_34, y=-log10(p12_34))) + 
  geom_point(alpha = 1, 
             shape = 16,
             size = 1, colour = '#00A087B2') +
  annotate("text", x = 1.6, y = 14, label = "189 sig. & logFC>0", 
           family= theme_get()$text[["family"]], )+
  annotate("text", x = -1.6, y = 14, label = "10 sig. & logFC<0",
           family= theme_get()$text[["family"]], )+
  xlab("log fold change") +
  ylab('-log10(FDR corrected p value)')+
  ggtitle('Expression associated with 1>=y>=2 vs y>=3')+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = c(seq(-3, 3, 1)),     
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(seq(0, 15, 5)),     
                     limits = c(0, 15))   +
  theme_Pub()


pdf('volcano p vs LFC 12_34 PBMC_sort global filtering .pdf', width = 3.5, height = 3.5)
plot(p1)
dev.off()


# sel = 

countdf[countdf$p12_34 == sort(countdf$p12_34)[7],]
peak_gene_pairs[5582,]




## 1 vs 2 

peak_gene_pairs[5571,]
# peak  gene
# 7778 chr6:144149912-144150816 STX11
peak_gene_pairs[3386,]
# peak  gene
# 3386 chr19:41876957-41877634 CD79A
peak_gene_pairs[1008,]
# peak gene
# 1008 chr11:10293701-10294466 SBF2
peak_gene_pairs[3395,]
# peak  gene
# 3395 chr19:43669865-43670691 PLAUR  -- inverse trend
peak_gene_pairs[2571,]
# peak gene
# 2571 chr17:7579041-7579956 CD68


vp = data.frame(pk = Xc_sub['chr6:150599271-150600181',], gene = Xr_sub['PLEKHG1',])

vp$pk[vp$pk >= 7] = 7
vp2 = vp[vp$pk %in% c(0,1,2,3,4,5,6,7),]
# vp2$pk[vp2$pk == 102] = 2
# vp2$pk[vp2$pk == 103] = 3
vp2$pk = as.character(vp2$pk)
vp2$pk[vp2$pk == '7'] = '>=7'
vp2$pk = factor(vp2$pk, levels = c('0', '1','2','3','4','5','6','>=7'))



p = ggplot(vp2, aes(x=pk, y=gene, fill = pk)) + 
  geom_boxplot(outlier.size=0.8) + 
  xlab("peak accessibility (Insertion)") +
  ylab('gene expression')+
  ggtitle('chr6:150599271-150600181 and PLEKHG1')+
  scale_y_continuous(breaks = c(seq(0, 40, 10)),     
                     limits = c(0, 40))   +
  theme(panel.background = element_blank(),legend.position= "none",
        plot.background = element_rect(colour = NA),
        plot.title = element_text(size=10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black", size = 0.2),
        axis.ticks = element_line(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        text = element_text(size = 10),
        axis.title = element_text(size=10, colour = 'black'),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  scale_fill_brewer(palette="Blues")

pdf('violing plot PLEKHG1 1t7.pdf',width = 3.5, height = 3.5)
plot(p)
dev.off()



p = ggplot(vp2, aes(x=pk, y=gene, fill = pk)) + 
  geom_boxplot(outlier.size=0.8) + 
  xlab("peak accessibility (Insertion)") +
  ylab('gene expression')+
  ggtitle('chr22:50529672-50530544 and TYMP')+
  scale_y_continuous(breaks = c(seq(0, 40, 10)),     
                     limits = c(0, 40))   +
  theme(panel.background = element_blank(),legend.position= "none",
        plot.background = element_rect(colour = NA),
        plot.title = element_text(size=10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black", size = 0.2),
        axis.ticks = element_line(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        text = element_text(size = 10),
        axis.title = element_text(size=10, colour = 'black'),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  scale_fill_brewer(palette="Blues")

pdf('violing plot TYMP 1t7.pdf',width = 3.5, height = 3.5)
plot(p)
dev.off()



p = ggplot(vp2, aes(x=pk, y=gene, fill = pk)) + 
  geom_boxplot(outlier.size=0.8) + 
  xlab("peak accessibility (Insertion)") +
  ylab('gene expression')+
  ggtitle('chr5:150412442-150413364 and CD74')+
  scale_y_continuous(breaks = c(seq(0, 40, 10)),     
                     limits = c(0, 40))   +
  theme(panel.background = element_blank(),legend.position= "none",
        plot.background = element_rect(colour = NA),
        plot.title = element_text(size=10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black", size = 0.2),
        axis.ticks = element_line(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        text = element_text(size = 10),
        axis.title = element_text(size=10, colour = 'black'),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  scale_fill_brewer(palette="Blues")

pdf('violing plot CD74 1t7.pdf',width = 3.5, height = 3.5)
plot(p)
dev.off()



p = ggplot(vp2, aes(x=pk, y=gene, fill = pk)) + 
  geom_boxplot(outlier.size=0.8) + 
  xlab("peak accessibility (Insertion)") +
  ylab('gene expression')+
  ggtitle('chr11:10293701-10294466 and SBF2')+
  scale_y_continuous(breaks = c(seq(0, 40, 10)),     
                     limits = c(0, 40))   +
  theme(panel.background = element_blank(),legend.position= "none",
        plot.background = element_rect(colour = NA),
        plot.title = element_text(size=10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black", size = 0.2),
        axis.ticks = element_line(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        text = element_text(size = 10),
        axis.title = element_text(size=10, colour = 'black'),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  scale_fill_brewer(palette="Blues")

pdf('violing plot SBF2 1t5.pdf',width = 3.5, height = 3.5)
plot(p)
dev.off()





count_sum <- matrix(nrow = dim(Xc_sub)[1], ncol = 15)
colnames(count_sum) <- c('mean_1', 'mean_2', 'mean_102', 'mean_2or102', 'mean_g103',
                         'nz_prop_1','nz_prop_2','nz_prop_102','nz_prop_2or102','nz_prop_g103',
                         'p1_2','p1_102','p1_2or102','p2_102','p102_g103')

gene_recorded <- vector(length = dim(Xc_sub)[1],mode = 'character')
peak_recorded <- rownames(Xc_sub)




## save Xc as full matrix
Xc_sub <- as.matrix(Xc_sub)
pk2_sub <- peak_sets

for(i in 1:dim(Xc_sub)[1]){
  atac_i = Xc_sub[i,]
  pk_i = pk2_sub[i]
  prom_i =  subsetByOverlaps(prom_sub, pk_i)
  gene_recorded[i] = prom_i$gene_name
  if( gene_recorded[i] %in% rownames(Xr_sub)){
    rna_i = Xr_sub[gene_recorded[i],]
    count_sum[i,] <- comp_summ_b(atac_i, rna_i)
  }else{
    next
  }
}

peak_gene_pairs = data.frame(peak = peak_recorded, gene = gene_recorded)

saveRDS(count_sum, 'count_sum_peak_gene_PBMC_10k_bi_allele_2or102_g103_sorted.rds')
saveRDS(peak_gene_pairs, 'peak_gene_pairs.rds')

count_sum = readRDS('count_sum_peak_gene_PBMC_10k_bi_allele_2or102_g103_sorted.rds')
peak_gene_pairs = readRDS('peak_gene_pairs.rds')

quantile(count_sum[,'p1_2'], na.rm = T)
hist(count_sum[,'p1_2'])
hist(count_sum[,'p1_102'])
hist(count_sum[,'p1_2or102'])
hist(count_sum[,'p102_g103'])
sum(count_sum[,'p1_2'] < 0.05, na.rm = T) 
sum(count_sum[,'p102_g103'] < 0.05, na.rm = T) 


cdf = as.data.frame(count_sum)
cdf$fc1_2 = log(cdf$mean_2 / cdf$mean_1)
cdf$fc1_102 = log(cdf$mean_102 / cdf$mean_1)
cdf$fc1_2or102 = log(cdf$mean_2or102 / cdf$mean_1)
cdf$fc2_102 = log(cdf$mean_102 / cdf$mean_2)
cdf$fc102_g103 = log(cdf$mean_g103 / cdf$mean_102)

hist(cdf[,'fc1_2'])
hist(cdf[,'fc1_102'])
hist(cdf[,'fc1_2or102'])
hist(cdf[,'fc102_g103'])



countdf = cbind(cdf, peak_gene_pairs)

## remove any row containing any NA
countdf = na.omit(countdf)
countdf$p1_2 = p.adjust(countdf$p1_2,method = 'fdr')
countdf$p1_102 = p.adjust(countdf$p1_102,method = 'fdr')
countdf$p2_102 = p.adjust(countdf$p2_102,method = 'fdr')
countdf$p102_g103 = p.adjust(countdf$p102_g103,method = 'fdr')
countdf$p1_2or102 = p.adjust(countdf$p1_2or102,method = 'fdr')



sum(countdf$p1_2or102 < 0.05 & countdf$fc1_2or102 > 0 )
sum(countdf$p1_2or102 < 0.05 & countdf$fc1_2or102 < 0 )

p1 <- ggplot(data=countdf, aes(x=fc1_2or102, y=-log10(p1_2or102))) + 
  geom_point(alpha = 1, 
             shape = 16,
             size = 1, colour = '#00A087B2') +
  annotate("text", x = 1.5, y = 10, label = "63 sig. & logFC>0", 
           family= theme_get()$text[["family"]], )+
  annotate("text", x = -1.5, y = 10, label = "5 sig. & logFC>0",
           family= theme_get()$text[["family"]], )+
  xlab("log fold change") +
  ylab('-log10(FDR corrected p value)')+
  ggtitle('Fragment-based counting, 1 vs 2')+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = c(seq(-3, 3, 1)),     
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(seq(0, 15, 5)),     
                     limits = c(0, 15))   +
  theme_Pub()


pdf('volcano FDR p vs LFC 1_2or102_10K_sorted rearragne small.pdf', width = 3.5, height = 3.5)
plot(p1)
dev.off()

countdf[countdf$p1_2or102 == min(countdf$p1_2or102, na.rm = T),]
countdf[countdf$fc1_2or102 == max(countdf$fc1_2or102, na.rm = T),]


## violin plot
vp = data.frame(pk = Xc_sub['chr19:41876957-41877634',], gene = Xr_sub['CD79A',])
vp2 = vp[vp$pk %in% c(0,1,2,102,3,103),]
vp2$pk[vp2$pk == 102] = 2
vp2$pk[vp2$pk == 103] = 3

vp2$pk = factor(vp2$pk, levels = c('0', '1','2','3'))


p = ggplot(vp2, aes(x=pk, y=gene,fill = pk)) + 
  geom_boxplot(outlier.size=0.8) + 
  xlab("peak accessibility") +
  ylab('gene expression')+
  ggtitle('chr19:41876957-41877634 and CD79A')+
  theme_Pub()+
  scale_fill_brewer(palette="Blues")

pdf('violing plot CD79A 102t2.pdf',width = 3.5, height = 4)
plot(p)
dev.off()

countdf[countdf$p1_2or102 == sort(countdf$p1_2or102)[2],]


vp = data.frame(pk = Xc_sub['chr6:144149912-144150816',], gene = Xr_sub['STX11',])
vp2 = vp[vp$pk %in% c(0,1,2,102,3,103),]
vp2$pk[vp2$pk == 102] = 2
vp2$pk[vp2$pk == 103] = 3

vp2$pk = factor(vp2$pk, levels = c('0', '1','2','3'))


p = ggplot(vp2, aes(x=pk, y=gene, fill = pk)) + 
  geom_boxplot(outlier.size=0.8) + 
  xlab("peak accessibility") +
  ylab('gene expression')+
  ggtitle('chr6:144149912-144150816 and STX11')+
  scale_y_continuous(breaks = c(seq(0, 40, 10)),     
                     limits = c(0, 40))   +
  theme(panel.background = element_blank(),legend.position= "none",
        plot.background = element_rect(colour = NA),
        plot.title = element_text(size=10,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black", size = 0.2),
        axis.ticks = element_line(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        text = element_text(size = 10),
        axis.title = element_text(size=10, colour = 'black'),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  scale_fill_brewer(palette="Blues")

pdf('violing plot STX11 102t2.pdf',width = 3.5, height = 4)
plot(p)
dev.off()

sum(countdf$p1_102 < 0.05 & countdf$fc1_102 > 0 )
sum(countdf$p1_102 < 0.05 & countdf$fc1_102 < 0 )




p1 <- ggplot(data=countdf, aes(x=fc1_102, y=-log10(p1_102))) + 
  geom_point(alpha = 1, 
             shape = 16,
             size = 1, colour = '#00A087B2') +
  annotate("text", x = 2.5, y = 8, label = "352 sig. & logFC>0", 
           family= theme_get()$text[["family"]], )+
  annotate("text", x = -2.5, y = 8, label = "53 sig. & logFC>0",
           family= theme_get()$text[["family"]], )+
  xlab("log fold change") +
  ylab('-log10(p value)')+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),     
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 10, 2)),     
                     limits = c(0, 10))   +
  theme_Pub()


pdf('volcano p vs LFC 1_102_10K_sorted.pdf', width = 4, height = 4)
plot(p1)
dev.off()


sum(countdf$p102_103 < 0.05 & countdf$fc102_103 > 0 )
sum(countdf$p102_103 < 0.05 & countdf$fc102_103 < 0 )


p1 <- ggplot(data=countdf, aes(x=fc102_103, y=-log10(p102_103))) + 
  geom_point(alpha = 1, 
             shape = 16,
             size = 1, colour = '#00A087B2') +
  annotate("text", x = 2.5, y = 8, label = "215 sig. & logFC>0", 
           family= theme_get()$text[["family"]], )+
  annotate("text", x = -2.5, y = 8, label = "99 sig. & logFC>0",
           family= theme_get()$text[["family"]], )+
  xlab("log fold change") +
  ylab('-log10(p value)')+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),     
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 10, 2)),     
                     limits = c(0, 10))   +
  theme_Pub()+
  scale_colour_aaas()


pdf('volcano p vs LFC 102_103_10K_sorted.pdf', width = 4, height = 4)
plot(p1)
dev.off()





p2 <- ggplot(data=countdf, aes(x=fc1_102, y=-log10(p1_102))) + 
  geom_point(alpha = 1, 
             shape = 16,
             size = 1,
             colour = '#4DBBD5B2') + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),     
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 15, 5)),     
                     limits = c(0, 15))   +
  theme_minimal()

plot(p2)


