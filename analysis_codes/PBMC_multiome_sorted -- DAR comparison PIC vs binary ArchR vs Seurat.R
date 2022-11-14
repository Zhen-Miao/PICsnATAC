

## load libraries
library('presto')
setwd('/home/zhenmiao/test_github')
## load a wrapper function for Seurat or ArchR approach
## this file below is also in the analysis_code folder
source('other_methods_for_differential cleaned.R')


## load necessary data
cell_type_labels <- readRDS('./pbmc2_CD4_CD14_subset_cluster_label.rds')
cell_type_labels <- as.character(cell_type_labels)
cell_type_set <- as.character(unique(cell_type_labels))

kn <- readRDS('./pbmc2_CD4_CD14_subset_ATAC_pic.rds')

## firstly, remove peaks with all zero
rskn = rowSums(kn)
quantile(rskn)
quantile(rskn,c(0.05,0.1,0.2))


## filter
kns = kn[rskn >= 20,]  # 130241 peaks * 4310 cells , 0.2491739 counts >=2
ctsub = cell_type_labels
kns@x[kns@x >= 7] <- 7

rs_kns_sum = readRDS('./row_sum_count_ge_2_BMMC_sorted_CD4_CD14.rds')
kns = kns[rs_kns_sum >= 10,]

knb = kns
knb@x = rep(1, length = length(knb@x))

dim(kns)
# [1] 61381  4310


mpos = readRDS('./mat1_pos_CD4.rds')
mneg = readRDS('./mat2_neg_CD14.rds')

p_archR <- archR_method_custom(data_matrix_pos_count = mpos,
                               data_matrix_neg_count = mneg,
                               k = 20,
                               bufferRatio = 1,
                               maxCells = 9000)
saveRDS(object = p_archR,file = 'p_archR_pic_CD4_CD14_all_cells.rds')

mpos = mpos
mpos@x <- rep(1, length = length(mpos@x))
mneg = mneg
mneg@x <- rep(1, length = length(mneg@x))

p_archR_b <- archR_method_custom(data_matrix_pos_count = mpos,
                               data_matrix_neg_count = mneg,
                               k = 20,
                               bufferRatio = 1,
                               maxCells = 9000)
# p_archR_b <- df$pval
saveRDS(object = p_archR_b,file = 'p_archR_pic_binary_CD4_CD14_all_cells.rds')


mpos = readRDS('./mat1_pos_CD4.rds')
mneg = readRDS('./mat2_neg_CD14.rds')
mpos = as.matrix(mpos)
mneg = as.matrix(mneg)
pfrag = c(Rfast::colsums(mpos), Rfast::colsums(mneg))
p_seurat = seurat_method2_subsample(mpos, mneg, pfrag)

saveRDS(object = p_seurat, file = 'p_seurat_pic_CD4_CD14_all_cells.rds')


mpos@x <- rep(1, length = length(mpos@x))
mneg@x <- rep(1, length = length(mneg@x))
mpos = as.matrix(mpos)
mneg = as.matrix(mneg)
pfrag = c(Rfast::colsums(mpos), Rfast::colsums(mneg))
p_seurat_b = seurat_method2_subsample(mpos, mneg, pfrag)

saveRDS(object = p_seurat_b, file = 'p_seurat_pic_binary_CD4_CD14_all_cells.rds')



## repeat 5 times, each time sample 1000 cells
set.seed(245)
sample_size = 500

mpos = readRDS('./mat1_pos_CD4.rds')
mneg = readRDS('./mat2_neg_CD14.rds')

n_np = dim(mpos)[2]
n_stroma = dim(mneg)[2]

for(s_round in c(1:5)){

  mpos = readRDS('./mat1_pos_CD4.rds')
  mneg = readRDS('./mat2_neg_CD14.rds')

  mpos@x = rep(1, length = length(mpos@x))
  mneg@x = rep(1, length = length(mneg@x))

  set.seed(245+s_round)
  np_sam = sample(x = n_np, size = sample_size,replace = F)
  st_sam = sample(x = n_stroma, size = sample_size,replace = F)

  mpos = mpos[,np_sam]
  mneg = mneg[,st_sam]

  p_archR <- archR_method_custom(data_matrix_pos_count = mpos,
                                 data_matrix_neg_count = mneg,
                                 k = 20,
                                 bufferRatio = 1,
                                 maxCells = 9000)
  saveRDS(object = p_archR,file =
            paste('p_archR_pic_binary_CD4_CD14_size',sample_size,'round',s_round,'.rds',sep = '_' ))
  mpos = as.matrix(mpos)
  mneg = as.matrix(mneg)
  pfrag = c(Rfast::colsums(mpos), Rfast::colsums(mneg))
  p_seurat = seurat_method2_subsample(mpos, mneg, pfrag)

  saveRDS(object = p_seurat, file =
          paste('p_seurat_pic_binary_CD4_CD14_size',sample_size,'round',s_round,'.rds',sep = '_' ))

}



