
## pbmc

# library("Seurat")
library(viridisLite);
library(ggplot2);

library("SingleCellExperiment")
library("monocle3")
library("cicero")
library("stringr")
library("FNN")
library("glasso")
library(dplyr)
library(Matrix)


cell_type_labels <- readRDS('./pbmc2_CD4_CD14_subset_cluster_label.rds')
cell_sel = names(cell_type_labels[cell_type_labels == 'CD14 Mono'])

all_cellinfo.8 = readRDS('./all_cellinfo.8_pbmc.rds')
all_peak_info.8 = readRDS('./all_peak_info.8_pbmc.rds')
all_mat.8 = readRDS('./all_mat.8_pbmc.rds')
umap_coords = readRDS('./umap_coords_pbmc.rds')

all_cellinfo.8 = all_cellinfo.8[cell_sel,]
umap_coords = umap_coords[cell_sel,]
all_mat.8 = all_mat.8[,cell_sel]

dim(all_mat.8)
# [1] 143823   2899

bc = read.table('./pbmc_10k_sorted_barcodes.tsv', sep = '\t')
bc = bc$V1
pk = read.table('./pbmc_10k_sorted_peaks.bed')
pk = pk$V1

inser_mat = readMM('./pbmc_10k_sorted_insertion_mat.mtx')
rownames(inser_mat) = pk
colnames(inser_mat) = bc


inser_mat.8 = inser_mat[rownames(all_mat.8), colnames(all_mat.8)]



input_cds_all.8 <-  suppressWarnings(new_cell_data_set(inser_mat.8,
                                                       cell_metadata = all_cellinfo.8,
                                                       gene_metadata = all_peak_info.8))

input_cds_all.8 <- monocle3::detect_genes(input_cds_all.8)


#Ensure there are no peaks included with zero reads
input_cds_all.8 <- input_cds_all.8[Matrix::rowSums(exprs(input_cds_all.8)) != 0,]


cicero_cds <- make_cicero_cds(input_cds_all.8, reduced_coordinates = umap_coords)



# input_cds.8 = input_cds_all.8

hg38 = data.frame(V1 = paste('chr', c(1:22, 'X', 'Y', 'M'), sep = ''),
                  V2 = c(248956422,242193529,198295559,190214555,181538259,
                         170805979,159345973,145138636,138394717,133797422,
                         135086622,133275309,114364328,107043718,101991189,
                         90338345,83257441,80373285,58617616,64444167,46709983,
                         50818468,156040895,	57227415,16569))

sample_genome <- subset(hg38, V1 == "chr1")
conns <- run_cicero(cicero_cds, sample_genome, sample_num = 100)
saveRDS(conns, 'conns_insertion_data_chr1_CD14.rds')

head(conns)


## binary
all_mat.8@x = rep(1, length = length(all_mat.8@x))

input_cds_all.8 <-  suppressWarnings(new_cell_data_set(all_mat.8,
                                                       cell_metadata = all_cellinfo.8,
                                                       gene_metadata = all_peak_info.8))

input_cds_all.8 <- monocle3::detect_genes(input_cds_all.8)


#Ensure there are no peaks included with zero reads
input_cds_all.8 <- input_cds_all.8[Matrix::rowSums(exprs(input_cds_all.8)) != 0,]


cicero_cds <- make_cicero_cds(input_cds_all.8, reduced_coordinates = umap_coords)


sample_genome <- subset(hg38, V1 == "chr1")
conns <- run_cicero(cicero_cds, sample_genome, sample_num = 100)
saveRDS(conns, 'conns_binary_data_chr1_CD14.rds')



