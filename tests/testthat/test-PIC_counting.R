test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("data.table is loadable", {
  library('PICsnATAC')

  meta.data <- read.csv('data/atac_singlecell_test_sub.csv', header = T )
  meta.data_filtered <- meta.data[meta.data$is__cell_barcode == 1,]
  cells <- meta.data_filtered$barcode

  peaks <- data.table::fread('data/atac_pbmc_5k_nextgem_peaks.bed',header = F)
  colnames(peaks) <- c('seqname', 'start', 'end')
  peak_sets = GenomicRanges::makeGRangesFromDataFrame(peaks)

  fragment_tsv_gz_file_location <- 'data/atac_fragment_test.tsv.gz'
  pic_mat <- PIC_counting(
    cells = cells,
    fragment_tsv_gz_file_location = fragment_tsv_gz_file_location,
    peak_sets = peak_sets)

  expect_equal(pic_mat@i, c(8, 8, 3, 1, 4, 6, 3, 7, 8, 4))
  expect_equal(pic_mat@x, c(1, 3, 1, 1, 1, 1, 1, 1, 1, 1))

})


