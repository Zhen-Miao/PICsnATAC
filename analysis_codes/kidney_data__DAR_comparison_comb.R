## test for type1 error simulation 

## load libraries
# library('SummarizedExperiment')
# library(aricode)
library('presto')
source('param_estimate_joint_p_q.R')
source('CCT_function.R')
source('differential_identification.R')
source('param_estimate_simple.R')
source('param-estimate_logit.R')
source('firth_implement_functions newton_try.R')
source('other_methods_for_differential.R')
source('param_estimate_doublet.R')

## load necessary data

cell_type_labels <- readRDS('kidcluster.rds')
cell_type_labels <- as.character(cell_type_labels)
samples <- readRDS('kidsample.rds')
cell_type_set <- as.character(unique(cell_type_labels))
peaks <- readRDS('peak_list_P0.rds')
cells <- readRDS('cells_sel.rds')

o28 = readRDS('peak by cell matrix insertion_kidney28.rds')
o29 = readRDS('peak by cell matrix insertion_kidney29.rds')

## now, match the cell types for the two groups
cells <- paste(cells, samples, sep = '_')
colnames(o28) <- paste(colnames(o28), '90028', sep = '_')
colnames(o29) <- paste(colnames(o29), '90029', sep = '_')

kn <- cbind(o28,o29)
kn <- kn[,cells]

## do DAR for NP and stroma2 clusters
kns = kn[,cell_type_labels %in% c('NP', 'stroma2')]
ctsub = cell_type_labels[cell_type_labels %in% c('NP', 'stroma2')]
kns@x[kns@x >= 14] <- 0 
kns <- Matrix::drop0(kns)

knb = kns
knb@x = rep(1, length = length(knb@x))
n_cell_acc <- Matrix::rowSums(knb)


o28f = readRDS('mtx_signac_P0_28.rds')
o29f = readRDS('mtx_signac_P0_29.rds')

## now, match the cell types for the two groups
# cells <- paste(cells, samples, sep = '_')
colnames(o28f) <- paste(colnames(o28f), '90028', sep = '_')
colnames(o29f) <- paste(colnames(o29f), '90029', sep = '_')

knf <- cbind(o28f,o29f)
knf <- knf[,cells]

## do DAR for NP and stroma2 clusters
knsf = knf[,cell_type_labels %in% c('NP', 'stroma2')]
ctsubf = cell_type_labels[cell_type_labels %in% c('NP', 'stroma2')]
knsf@x[knsf@x >= 7] <- 0 
knsf <- Matrix::drop0(knsf)

knbf = knsf
knbf@x = rep(1, length = length(knbf@x))
n_cell_accf <- Matrix::rowSums(knbf)

ss = n_cell_acc >= 20
jj = n_cell_accf >= 20
sel_in = ss[jj]
sel_fr = jj[ss]


# should use the same filtration as

kns <- kns[ss & jj,]
knsf <- knsf[ss & jj,]

dim(kns)
# [1] 199464   4488
mpos = as.matrix(kns[,ctsub == 'NP'])
mneg = as.matrix(kns[,ctsub == 'stroma2'])
pfrag = c(Rfast::colsums(mpos), Rfast::colsums(mneg))
p_seurat = seurat_method2_subsample(mpos, mneg, pfrag)

saveRDS(p_seurat, 'p_seurat_insertion_NP_stroma2 true_fragment_mtx.rds')

p_archR = archR_method(data_matrix_pos_count = mpos,data_matrix_neg_count = mneg)
saveRDS(p_archR,'p_archR_insertion_NP_stroma2 true_fragment_mtx.rds' )

table(cell_type_labels)


mpos = as.matrix(knsf[,ctsubf == 'NP'])
mneg = as.matrix(knsf[,ctsubf == 'stroma2'])
pfrag = c(Rfast::colsums(mpos), Rfast::colsums(mneg))

p_seurat = seurat_method2_subsample(mpos, mneg, pfrag)
saveRDS(p_seurat, 'p_seurat_fragment_NP_stroma2 true_fragment_mtx.rds')

p_archR = archR_method(data_matrix_pos_count = mpos,data_matrix_neg_count = mneg)
saveRDS(p_archR,'p_archR_fragment_NP_stroma2 true_fragment_mtx.rds' )

table(cell_type_labels)




parchri = readRDS('p_archR_insertion_NP_stroma2.rds')
parchrf = readRDS('p_archR_fragment_NP_stroma2 true_fragment_mtx.rds')
# parchri = parchri[sel_fr]
# parchrf = parchrf[sel_in]

parch = data.frame(insertion = parchri, fragment = parchrf)
parchfdr = parch
parchfdr$insertion = p.adjust(parchfdr$insertion, method = 'fdr')
parchfdr$fragment = p.adjust(parchfdr$fragment, method = 'fdr')

sum(parchfdr$insertion < 0.05 & parchfdr$fragment < 0.05) # 141717
sum(parchfdr$insertion < 0.05 & parchfdr$fragment > 0.05) # 1043
sum(parchfdr$insertion > 0.05 & parchfdr$fragment < 0.05) # 1270


pseui = readRDS('p_seurat_insertion_NP_stroma2.rds')
pseuf = readRDS('p_seurat_fragment_NP_stroma2 true_fragment_mtx.rds')
pseui = pseui[sel_fr]
# pseuf = pseuf[sel_in]

pse = data.frame(insertion = pseui, fragment = pseuf)
psefdr = pse
psefdr$insertion = p.adjust(psefdr$insertion, method = 'fdr')
psefdr$fragment = p.adjust(psefdr$fragment, method = 'fdr')

sum(psefdr$insertion < 0.05 & psefdr$fragment < 0.05) # 140118
sum(psefdr$insertion < 0.05 & psefdr$fragment > 0.05) # 2689
sum(psefdr$insertion > 0.05 & psefdr$fragment < 0.05) # 4290

## 195513 total peaks tested 


r_by_ct_est <-  get_r_by_ct_mat_pq(
  cell_type_set = cell_type_set,
  r_by_c = kn,
  cell_type_labels = cell_type_labels,
  n_features_per_cell = dim(kn)[1],
  p_acc = 0.0005,
  q_acc = 0.0005,
  n_max_iter = 500
)

saveRDS(r_by_ct_est, 'r_by_ct_est_P0kidney.rds')

