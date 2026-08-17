test_that("PIC_counting reproduces the reference in-memory result", {
  meta_data <- utils::read.csv(
    testthat::test_path("data", "atac_singlecell_test_sub.csv")
  )
  cells <- meta_data$barcode[meta_data$is__cell_barcode == 1]

  peaks <- data.table::fread(
    testthat::test_path("data", "atac_pbmc_5k_nextgem_peaks.bed"),
    header = FALSE
  )
  colnames(peaks) <- c("seqname", "start", "end")
  peak_sets <- GenomicRanges::makeGRangesFromDataFrame(peaks)

  pic_mat <- PIC_counting(
    cells = cells,
    fragment_tsv_gz_file_location = testthat::test_path(
      "data", "atac_fragment_test.tsv.gz"
    ),
    peak_sets = peak_sets,
    verbose = FALSE
  )

  expect_s4_class(pic_mat, "sparseMatrix")
  expect_equal(pic_mat@i, c(8, 8, 3, 1, 4, 6, 3, 7, 8, 4))
  expect_equal(pic_mat@x, c(1, 3, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_identical(colnames(pic_mat), cells)
})


test_that("streamed and in-memory counting agree in arbitrary peak order", {
  skip_if_not_installed("Rsamtools")

  fragment_file <- tempfile(fileext = ".tsv")
  bgzip_file <- tempfile(fileext = ".tsv.gz")
  writeLines(c(
    "chr1\t99\t110\tcell-a\t1",
    "chr1\t299\t310\tcell-b\t1",
    "chr2\t199\t220\tcell-a\t1"
  ), fragment_file)
  Rsamtools::bgzip(fragment_file, dest = bgzip_file, overwrite = TRUE)
  Rsamtools::indexTabix(bgzip_file, format = "bed")
  on.exit(unlink(c(fragment_file, bgzip_file, paste0(bgzip_file, ".tbi"))))

  peaks <- data.frame(
    seqname = c("chr2", "chr-missing", "chr1", "chr1"),
    start = c(200, 1, 100, 300),
    end = c(210, 20, 105, 305)
  )
  cells <- c("cell-a", "cell-b", "cell-with-no-fragments")

  full <- suppressWarnings(PIC_counting(
    cells = cells,
    fragment_tsv_gz_file_location = bgzip_file,
    peak_sets = peaks,
    load_full = TRUE,
    verbose = FALSE
  ))
  streamed <- PIC_counting(
    cells = cells,
    fragment_tsv_gz_file_location = bgzip_file,
    peak_sets = peaks,
    load_full = FALSE,
    verbose = FALSE
  )

  expect_equal(as.matrix(streamed), as.matrix(full))
  expect_identical(rownames(streamed), rownames(full))
  expect_identical(colnames(streamed), cells)
  expect_true(all(streamed[2, ] == 0))
})


test_that("PIC_counting validates structural inputs", {
  peaks <- data.frame(seqname = "chr1", start = 1, end = 10)
  fragment_file <- testthat::test_path("data", "atac_fragment_test.tsv.gz")

  expect_error(
    PIC_counting(c("a", "a"), fragment_file, peaks),
    "duplicate"
  )
  expect_error(
    PIC_counting("a", fragment_file, peaks, extend_size = 1.5),
    "non-negative integer"
  )
  expect_error(data_frame_to_GRanges(data.frame(x = 1, y = 2)),
    "at least three"
  )
})
