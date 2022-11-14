## test for type1 error simulation

library(ggsci)
library('Matrix')
library('dplyr')

kns = readRDS('./pic_mat_CD4_CD14_filterred.rds')
capt_rate = readRDS('./cap_rate_CD4_CD14_filterred.rds')
ctsub = readRDS('./ctsub_CD4_CD14_filterred.rds')
p = readRDS('./p_val_CD4_CD14_MLE_all.rds')
parch = readRDS('./p_archR_pic_CD4_CD14_all_filterred.rds')
pseu = readRDS('./p_seurat_pic_CD4_CD14_all_filterred.rds')


### define cutoff by the 5th rank percentile of all p values
pic_cut = readRDS('./pic_cut.rds')
pseu_cut_n = readRDS('./pseu_cut.rds')
parch_cut_n = readRDS('./parch_cut.rds')

sum(p.adjust(p, 'fdr') < pic_cut)
sum(p.adjust(pseu,'fdr') < pseu_cut_n)
sum(p.adjust(parch,'fdr') < parch_cut_n)


## downsampled
psub = readRDS('./p_val_CD4_CD14_all_filterred_sub_500_MLE_all.rds')
pseu_sub = readRDS('./p_seurat_pic_CD4_CD14_all_cells_filterred_sub_500.rds')
parch_sub = readRDS('./p_archR_pic_CD4_CD14_all_cells_filterred_sub_500.rds')

## unbalanced
psub_unbal = readRDS('p_val_CD4_CD14_all_filterred_sub_500_MLE_all_unbal_new_more_acc.rds')
pseu_unbal = readRDS('p_seurat_pic_CD4_CD14_all_cells_filterred_sub_500_unbal_new.rds')
parch_unbal = readRDS('p_archR_pic_CD4_CD14_all_cells_filterred_sub_500_unbal_new.rds')


dfb = data.frame(
  methods = rep(c('PIC model', 'Seurat', 'ArchR'), times = 3),
  settings = rep(c('full data', 'down_500', 'down_500_unbal'), each = 3),
  power = c(
    sum(p.adjust(p,method = 'fdr') < pic_cut),
    sum(p.adjust(pseu,method = 'fdr') < pseu_cut_n),
    sum(p.adjust(parch,method = 'fdr') < parch_cut_n),

    sum(p.adjust(psub,method = 'fdr') < pic_cut),
    sum(p.adjust(pseu_sub,method = 'fdr') < pseu_cut_n),
    sum(p.adjust(parch_sub,method = 'fdr') < parch_cut_n),

    sum(p.adjust(psub_unbal,method = 'fdr') < pic_cut),
    sum(p.adjust(pseu_unbal,method = 'fdr') < pseu_cut_n),
    sum(p.adjust(parch_unbal,method = 'fdr') < parch_cut_n)
  )
)


dfb$methods = factor(x = dfb$methods, levels = c('PIC model', 'Seurat', 'ArchR'))
dfb$settings = factor(x = dfb$settings, levels = c('full data', 'down_500', 'down_500_unbal'))


## missing corrected row mean
capturing_rates = capt_rate
kns_mat1 = kns[,ctsub == 'CD4 Naive']
kns_mat2 = kns[,ctsub == 'CD14 Mono']

capturing_rates = ceiling(capturing_rates*10)/10

rm1 = rowSums(kns_mat1 %*% diag(1 / capturing_rates[ctsub == 'CD4 Naive'])) / sum(ctsub == 'CD4 Naive')
rm2 = rowSums(kns_mat2 %*% diag(1 / capturing_rates[ctsub == 'CD14 Mono'])) / sum(ctsub == 'CD14 Mono')

rm1_s = rowSums(kns_mat1 ) / sum(ctsub == 'CD4 Naive')
rm2_s = rowSums(kns_mat2 ) / sum(ctsub == 'CD14 Mono')

log_fc = log(rm2_s / (rm1_s+0.0001))
log_fc_c = log(rm2 / (rm1+0.0001))


## redo the analysis, using the true union set as standard -- optionally, we can also filter by log fold change
## filtering by log fold change does not affect the results
union_true = (p.adjust(p,method = 'fdr') < pic_cut |
  p.adjust(pseu,method = 'fdr') < pseu_cut_n |
  p.adjust(parch,method = 'fdr') < parch_cut_n) & abs(log_fc_c) > 0.1

union_true = (p.adjust(p,method = 'fdr') < pic_cut |
                p.adjust(pseu,method = 'fdr') < pseu_cut_n |
                p.adjust(parch,method = 'fdr') < parch_cut_n)


dfb2 = data.frame(
  methods = rep(c('PIC model', 'Seurat', 'ArchR'), times = 3),
  settings = rep(c('full data', '500', '500_unbal'), each = 3),
  power = c(
    sum(p.adjust(p,method = 'fdr') < pic_cut) / sum(union_true),
    sum(p.adjust(pseu,method = 'fdr') < pseu_cut_n)/ sum(union_true),
    sum(p.adjust(parch,method = 'fdr') < parch_cut_n)/ sum(union_true),

    sum(p.adjust(psub,method = 'fdr') < pic_cut & union_true)/ sum(union_true),
    sum(p.adjust(pseu_sub,method = 'fdr') < pseu_cut_n & union_true)/ sum(union_true),
    sum(p.adjust(parch_sub,method = 'fdr') < parch_cut_n & union_true)/ sum(union_true),

    sum(p.adjust(psub_unbal,method = 'fdr') < pic_cut& union_true)/ sum(union_true),
    sum(p.adjust(pseu_unbal,method = 'fdr') < pseu_cut_n& union_true)/ sum(union_true),
    sum(p.adjust(parch_unbal,method = 'fdr') < parch_cut_n& union_true)/ sum(union_true)
  )
)

dfb2$methods = factor(x = dfb2$methods, levels = c('PIC model', 'Seurat', 'ArchR'))
dfb2$settings = factor(x = dfb2$settings, levels = c('full data', '500', '500_unbal'))





df = data.frame(ppic = p.adjust(p, 'fdr'),
                pseu = p.adjust(pseu, 'fdr'),
                parch = p.adjust(parch,'fdr'),
                rowmean1_c = rm1,rowmean2_c = rm2,
                rowmean1_s = rm1_s, rowmean2_s = rm2_s
)

## power vs row mean analysis

local_power = matrix(nrow = 25, ncol = 3)
# colnames(local_power) = c('PIC','Seurat', 'ArchR', 'PIC_500','Seurat_500','ArchR_500')
colnames(local_power) = c('PIC','Seurat', 'ArchR')
mean_rm_s = (rm1_s + rm2_s) / 2
qtls = quantile(mean_rm_s, probs = seq(0,1,0.04))
for(i in c(1:25)){
  pk_sel = which(mean_rm_s >= qtls[i] & mean_rm_s < qtls[i+1])
  n_true_pos = sum(union_true[pk_sel])
  local_power[i,'PIC'] = sum(df$ppic[pk_sel] < pic_cut ) / n_true_pos
  local_power[i,'Seurat'] = sum(df$pseu[pk_sel] < pseu_cut_n ) / n_true_pos
  local_power[i,'ArchR'] = sum(df$parch[pk_sel] < parch_cut_n ) / n_true_pos
}

dfpower = data.frame(
  local_power = as.vector(local_power),
  methods = rep(c('PIC model', 'Seurat', 'ArchR'), each = 25),
  percentiles = 1:25 * 4
)
dfpower$methods = factor(dfpower$methods, levels = c('PIC model', 'Seurat', 'ArchR'))



