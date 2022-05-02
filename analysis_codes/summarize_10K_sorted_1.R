
library(viridisLite);
library(ggplot2);
library(rtracklayer)
library(dplyr)
library('GenomicRanges')
library('Matrix')
library("ggsci")
library("gridExtra")


cells = readRDS('pbmc_granulocyte_sorted_10k_cells_wo_doublets.rds')
peak_sets = readRDS('peak_subset_in_gene_TSS.rds')
n_features = length(peak_sets)
pkdf = as.data.frame(peak_sets)
fname = paste(pkdf$seqnames, ':', pkdf$start, '-', pkdf$end, sep = '')

peak_ind = 1:n_features 

out_mat = readRDS('out_matrix_promoters_PBMC_10k_sort __ simple insertion.rds')
rownames(out_mat) = fname

inputdata.10x <- Seurat::Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
# atac_counts <- inputdata.10x$Peaks
rm(inputdata.10x)

## filter RNA based on selected cells 

rna_counts <- rna_counts[,cells]
rna_counts <- rna_counts[rowSums(rna_counts != 0) >= 5,]
# atac_counts <- atac_counts[,cells]

Xr <- apply(rna_counts, 2, function(x) (x/sum(x))*10000)
colSums(Xr[,1:5])


nexp = rowSums(Xr != 0)
quantile(nexp)
# 0%   25%   50%   75%  100%
# 5    25   174  1062 11214
Xr <- Xr[nexp >= 100,]   ## 13633 11234



## promoter overlap 
seqlevelsStyle(peak_sets) <- "NCBI"


gene.ranges = readRDS('gene_region_longest_HS_v86_signac.rds')
# tss <- resize(gene.ranges, width = 1, fix = 'start')
prom <- promoters(gene.ranges, upstream=100, downstream=100,use.names=TRUE)


prom_sub = subsetByOverlaps(prom, peak_sets)
gene_sel = prom_sub$gene_name


gene_sel_sec = intersect(rownames(Xr), gene_sel)  ## 5762


Xr_sub = Xr[gene_sel_sec,]
Xc_sub = out_mat[fname,]

rm(out_mat, Xr)


comp_summ_b <- function(
  atac_vec,
  rna_vec
){
  
  ret <- vector(length = 11)
  names(ret) <- c('mean_1', 'mean_2', 'mean_12', 'mean_34', 
                  'nz_prop_1','nz_prop_2','nz_prop_12','nz_prop_34','nz_prop_g0',
                  'p1_2','p12_34')
  
  ## firstly, require globally have 15% non-zero counts
  isg0 = (atac_vec >= 1)
  ret['nz_prop_g0'] = mean(rna_vec[isg0] != 0)
  
  is1 = (atac_vec == 1)
  is2 = (atac_vec == 2)
  is12 = (is1 | is2)
  is34 = (atac_vec >= 3 )
  
  ret['mean_1'] = mean(rna_vec[is1])
  ret['mean_2'] = mean(rna_vec[is2])
  ret['mean_12'] = mean(rna_vec[is12])
  ret['mean_34'] = mean(rna_vec[is34])
  
  ret['nz_prop_1'] = mean(rna_vec[is1] != 0)
  ret['nz_prop_2'] = mean(rna_vec[is2] != 0)
  ret['nz_prop_12'] = mean(rna_vec[is12] != 0)
  ret['nz_prop_34'] = mean(rna_vec[is34] != 0)
  
  ret['p1_2'] = wilcox.test(x = rna_vec[is1], y = rna_vec[is2], alternative = 'two.sided' )$p.value
  ret['p12_34'] = wilcox.test(x = rna_vec[is12], y = rna_vec[is34], alternative = 'two.sided' )$p.value
  return(ret)
}



count_sum <- matrix(nrow = dim(Xc_sub)[1], ncol = 11)
colnames(count_sum) <- c('mean_1', 'mean_2', 'mean_12', 'mean_34', 
                         'nz_prop_1','nz_prop_2','nz_prop_12','nz_prop_34','nz_prop_g0',
                         'p1_2','p12_34')

gene_recorded <- vector(length = dim(Xc_sub)[1],mode = 'character')
peak_recorded <- rownames(Xc_sub)

pk2_sub = peak_sets

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


saveRDS(count_sum, 'count_sum_peak_gene_PBMC_sorted_1234 g3 global_filtering .rds')
saveRDS(peak_gene_pairs, 'peak_gene_pairs_PBMC_sorted_1234 g3 global_filtering  .rds')
# saveRDS(count_sum, 'count_sum_peak_gene_PBMC_sorted_1234 g3 .rds')
# saveRDS(peak_gene_pairs, 'peak_gene_pairs_PBMC_sorted_1234 g3 .rds')
# 
# saveRDS(Xc_sub,'Xc_sub_1234PBMC_sorted.rds')
# saveRDS(Xr_sub,'Xr_sub_1234PBMC_sorted.rds')
