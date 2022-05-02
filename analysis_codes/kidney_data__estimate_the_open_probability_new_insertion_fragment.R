## test for type1 error simulation 

## load libraries
# library('SummarizedExperiment')
# library(aricode)
source('param_estimate_joint_p_q.R')
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

## --------------------- insertion --------------------- 

o28 = readRDS('peak by cell matrix insertion_kidney28.rds')
o29 = readRDS('peak by cell matrix insertion_kidney29.rds')
## now, match the cell types for the two groups
cells <- paste(cells, samples, sep = '_')
colnames(o28) <- paste(colnames(o28), '90028', sep = '_')
colnames(o29) <- paste(colnames(o29), '90029', sep = '_')
kn <- cbind(o28,o29)
kn <- kn[,cells]

## --------------------- fragment --------------------- 

o28f = readRDS('mtx_signac_P0_28.rds')
o29f = readRDS('mtx_signac_P0_29.rds')
colnames(o28f) <- paste(colnames(o28f), '90028', sep = '_')
colnames(o29f) <- paste(colnames(o29f), '90029', sep = '_')
knf <- cbind(o28f,o29f)
knf <- knf[,cells]

## --------------------- check the correspondence --------------------- 
## make sure the colnames match
all(grep(':','-', rownames(kn)) == rownames(knf)) ## T

ji = findInterval(seq(kn@x)-1,kn@p[-1])+1
jf = findInterval(seq(knf@x)-1,knf@p[-1])+1

dfi = data.frame(i = kn@i, j = ji, x1 = kn@x)
dff = data.frame(i = knf@i, j = jf, x2 = knf@x)

# df3 = merge(df1,df2, by = c('i', 'j'), all.x = F, all.y = F)
# df4 = inner_join(dfi, dff, by = c('i','j'))
# colnames(df4) = c('i','j', 'count_insertion', 'count_fragment')

df4 = full_join(dfi, dff,by = c('i','j'))
colnames(df4) = c('i','j', 'count_insertion', 'count_fragment')

sum(df4$count_insertion == 1 & is.na(df4$count_fragment))
sum(df4$count_fragment == 1 & is.na(df4$count_insertion))
df4$count_insertion[is.na(df4$count_insertion)] <- 0
df4$count_fragment[is.na(df4$count_fragment)] <- 0

sum(df4$count_fragment == 1 & df4$count_insertion == 1)
sum(df4$count_fragment == 1 & df4$count_insertion == 2)

sum(df4$count_fragment == 2 & df4$count_insertion == 0)
sum(df4$count_fragment == 2 & df4$count_insertion == 1)
sum(df4$count_fragment == 2 & df4$count_insertion == 2)


## --------------------- check open probability --------------------- 

## we should use insertion based counting, then, binarize it 

## remove peaks with too many cells getting >= 7
knn = knf
knn@x[knn@x <= 7 ] = 0
knn = drop0(knn)
knn@x = rep(1, length(knn@x))

n7 = rowSums(knn)
sel_n7 <- n7 <= 10  ## 37 peaks removed


kn@x[kn@x > 14] <- 0
kn <- drop0(kn)

knb = kn
knb@x = rep(1, length = length(knb@x))

n_cell_acc <- Matrix::rowSums(knb)
# 0%  25%  50%  75% 100%
# 7   46   95  224 6110
knb <- knb[n_cell_acc >= 40,]
sel_nc40 <- n_cell_acc >= 40


r_by_ct_est <-  get_r_by_ct_mat_pq(
  cell_type_set = cell_type_set,
  r_by_c = knb,
  cell_type_labels = cell_type_labels,
  n_features_per_cell = dim(knb)[1],
  p_acc = 0.0005,
  q_acc = 0.0005,
  n_max_iter = 500
)

saveRDS(r_by_ct_est, 'r_by_ct_est_P0kidney_new_insertion_counting.rds')

r_by_ct_est = readRDS('r_by_ct_est_P0kidney_new_insertion_counting.rds')


## ---------------------- get f1i1 and f1i2 ----------------------------

## filter peaks
r_by_ct_est$p_by_t_new = r_by_ct_est$p_by_t_new[sel_n7[sel_nc40],] 
kn = kn[sel_n7 & sel_nc40,]
knf = knf[sel_n7 & sel_nc40,]

## based on the insertion peak, check its correspondence 
ji = findInterval(seq(kn@x)-1,kn@p[-1])
jf = findInterval(seq(knf@x)-1,knf@p[-1])

dfi = data.frame(i = kn@i, j = ji, x1 = kn@x)
dff = data.frame(i = knf@i, j = jf, x2 = knf@x)

dfis = left_join(dfi,dff,by = c('i','j'))
sum(dfis$x1 == 1 & dfis$x2 == 1) # 12696631
sum(dfis$x1 == 2 & dfis$x2 == 1) # 40335039

dfis$x1[dfis$x1 == 1 & dfis$x2 == 1] <-  -1
dfis$x1[dfis$x1 == 2 & dfis$x2 == 1] <-  -2


## ---------------------- get f1i0 / f1i0 + f1i1 + f1i2 ----------------------------

## based on the insertion peak, check its correspondence 
ji = findInterval(seq(kn@x)-1,kn@p[-1])
jf = findInterval(seq(knf@x)-1,knf@p[-1])

dfi = data.frame(i = kn@i, j = ji, x1 = kn@x)
dff = data.frame(i = knf@i, j = jf, x2 = knf@x)

dfis = right_join(dfi,dff,by = c('i','j'))
dfis$x1[is.na(dfis$x1)] <- 0
sum(dfis$x1 == 1 & dfis$x2 == 1) # 12696631
sum(dfis$x1 == 2 & dfis$x2 == 1) # 40335039
sum(dfis$x1 == 0 & dfis$x2 == 1) # 269937

dfis$x2[dfis$x1 == 1 & dfis$x2 == 1] <-  -1
dfis$x2[dfis$x1 == 2 & dfis$x2 == 1] <-  -1
dfis$x2[dfis$x1 == 0 & dfis$x2 == 1] <-  -2




### ------------------- get the proportion of f1i2 and open prob correlation
pmatpt = sparseMatrix(i = dfis$i[dfis$x1<0], j = dfis$j[dfis$x1<0] , 
                      x = dfis$x1[dfis$x1<0],dims =dim(kn), index1=FALSE ) ## f1i2/f1i1+f1i2

pmatpt = sparseMatrix(i = dfis$i[dfis$x2<0], j = dfis$j[dfis$x2<0] , 
                      x = dfis$x2[dfis$x2<0],dims =dim(kn), index1=FALSE ) ## f1i2/f1i1+f1i2

pmatpt@x = -1 * pmatpt@x

## knf
pmatpt <- kn
pmatptbin <- pmatpt
pmatptbin@x <- rep(1, times = length(pmatptbin@x))

pmatpt@x[pmatpt@x == 2] = 1
pmatpt@x[pmatpt@x >= 3] = 2

ctypes <- cell_type_set 
proportion_counts_more_than_1 <- matrix(nrow = dim(pmatpt)[1], ncol = length(cell_type_set))
colnames(proportion_counts_more_than_1) <- cell_type_set

n_cell_acc_prop <- proportion_counts_more_than_1

## only look at 2 vs 1 
for(jj in ctypes){
  pmatpt1 = pmatpt[,cell_type_labels == jj]
  pmatptbin1 = pmatptbin[,cell_type_labels == jj]
  
  n_cell_acc <- Matrix::rowSums(pmatptbin1)
  n_counts_more_than_1 <- Matrix::rowSums(pmatpt1) - n_cell_acc
  proportion_counts_more_than_1[,jj] <- n_counts_more_than_1 / n_cell_acc
  n_cell_acc_prop[,jj] <- r_by_ct_est$p_by_t_new[,jj]
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



## for 3 vs 4
quantile(corr_across_features,na.rm = T)
# 0%        25%        50%        75%       100%
# -0.9642857  0.3571429  0.6123724  0.7783118  1.0000000
sum(corr_across_features < 0 , na.rm = T)  # 12948
sum(corr_across_features > 0 , na.rm = T) ## 171493
quantile(corr_p_across_features,na.rm = T)

sum(corr_p_across_features < 0.05, na.rm= T )  ## 57314
sum(corr_p_across_features < 0.05 & corr_across_features > 0, na.rm= T )  ## 57175

sel_p = corr_p_across_features[!is.na(corr_p_across_features) & corr_across_features > 0 ]
adj = p.adjust(sel_p)




df = data.frame(corr = corr_across_features)

p = ggplot(df, aes(x = corr)) + 
  geom_histogram(bins = 20, color='black', fill="#8491b4") + 
  xlim(-1,1)+
  # scale_x_discrete(drop = FALSE, limits = c(-1,1))+
  ggtitle("Correlation ( p, freq(y>=3|y>0) )") +
  xlab("Spearman Correlation") +
  ylab('Frequency')+
  theme_Pub()
plot(p)


pdf('frequency of correlation ig3 over ig0_kidney.pdf', width = 3.5, height = 3.5)
plot(p)
dev.off()






