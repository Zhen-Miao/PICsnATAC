

library('Matrix')


# get data

## load necessary data
cell_type_labels <- readRDS('./pbmc2_CD4_CD14_subset_cluster_label.rds')
cell_type_labels <- as.character(cell_type_labels)
kn <- readRDS('./pbmc2_CD4_CD14_subset_ATAC_pic.rds')

## set some values
n_cells_acc_filtering = 20
n_counts_g2_filtering = 10
max_pic_count_reasonable = 7

## filter
rskn = rowSums(kn)
pic_mat = kn[rskn >= 2,]  ## 143592   4310

# save some values
n_pks <- dim(pic_mat)[1]
ct_uniq <- unique(cell_type_labels)


## read saved capturing rate
r_by_ct_est = readRDS('./r_by_ct_est_CD4_CD14.rds')

capturing_rate_estimated = r_by_ct_est$q_vec_new
pic_b <- pic_mat
pic_b@x <- rep(1, length = length(pic_b@x))
rsb = rowSums(pic_b)
pic_f <- pic_mat[rsb >= n_cells_acc_filtering,]
pic_f@x[pic_f@x > max_pic_count_reasonable] <- max_pic_count_reasonable # 128866   4310

kns_sum <- pic_f
kns_sum@x[kns_sum@x == 1] <- 0
kns_sum <- drop0(kns_sum)
kns_sum@x <- rep(1, length = length(kns_sum@x))
rs_kns_sum <- rowSums(kns_sum)
pic_f = pic_f[rs_kns_sum >= n_counts_g2_filtering,] ## 61381  4310



rm(list = c('pic_mat','pic_b','kn'))


p_valall <- DAR_by_LRT(pic_mat = pic_f,
                     capturing_rates = capturing_rate_estimated,
                     cell_type_labels = cell_type_labels)

saveRDS(p_valall, 'p_val_CD4_CD14_bi_allele_all.rds')


## subset to 500 cells
pic_f1 = pic_f[,cell_type_labels == "CD4 Naive"]
pic_f2 = pic_f[,cell_type_labels == "CD14 Mono"]

capturing_rate_estimated_1 = capturing_rate_estimated[cell_type_labels == "CD4 Naive"]
capturing_rate_estimated_2 = capturing_rate_estimated[cell_type_labels == "CD14 Mono"]
set.seed(2655)
sf1 = sample(1:dim(pic_f1)[2], size = 500,replace = F)
sf2 = sample(1:dim(pic_f2)[2], size = 500,replace = F)


pic_f1 = pic_f1[,sf1]
pic_f2 = pic_f2[,sf1]

capturing_rate_estimated_1 =capturing_rate_estimated_1[sf1]
capturing_rate_estimated_2 =capturing_rate_estimated_2[sf2]

pic_fs = cbind(pic_f1,pic_f2)
cre_s = c(capturing_rate_estimated_1,capturing_rate_estimated_2)
cell_type_labels_s = rep(c("CD4 Naive","CD14 Mono"), each = 500)

rm(pic_f1,pic_f2)

p_valall500 <- DAR_by_LRT(pic_mat = pic_fs,
                        capturing_rates = cre_s,
                        cell_type_labels = cell_type_labels_s)

saveRDS(p_valall500, 'p_val_CD4_CD14_bi_allele_sub_500.rds')



## for random labeling to control type 1 error
set.seed(443 + 146)
cell_type_labels = sample(cell_type_labels,replace = F)


### binarize
pic_b <- pic_f
pic_b@x <- rep(1, length = length(pic_b@x))

### calculate p and q
r_by_ct_est <-  get_r_by_ct_mat_pq(
  cell_type_set = ct_uniq,
  r_by_c = pic_b,
  cell_type_labels = cell_type_labels,
  n_features_per_cell = dim(pic_b)[1]
)
pk_sample = sample(1:dim(pic_f)[1],size = 10000,replace = F)

p_val_permu2 <- DAR_by_LRT(pic_mat = pic_f[pk_sample,],
                     capturing_rates = r_by_ct_est$q_vec_new,
                     cell_type_labels = cell_type_labels)

saveRDS(p_val_permu2,'p_val_permu_r1.rds')


set.seed(27)
pk_sample = sample(1:dim(pic_f)[1],size = 10000,replace = F)

p_val_permu2 <- .DAR_by_LRT(pic_mat = pic_f[pk_sample,],
                            capturing_rates = r_by_ct_est$q_vec_new,
                            cell_type_labels = cell_type_labels)

saveRDS(p_val_permu2,'p_val_permu_r2.rds')

## repeat the above codes several times

## for gene activity score correlation analysis, do this
kn@x[kn@x > 7] <- 7 # 128866   4310
capturing_rates = r_by_ct_est$q_vec_new
pic_mat = kn

lamb_full <- obs_to_insertion_ME(pic_mat = pic_mat,
                                 capturing_rates = capturing_rates,
                                 cell_type_labels = cell_type_labels
)
## remove those not accessible
cd4 = pic_mat[,cell_type_labels == "CD4 Naive"]
cd4sum = rowSums(cd4)

cd14 = pic_mat[,cell_type_labels == "CD14 Mono"]
cd14sum = rowSums(cd14)


