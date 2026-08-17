test_that("DAR count preparation makes the artifact heuristic explicit", {
  dense <- matrix(c(0, 5, 6, 9, 10, 11, 15), nrow = 1)
  expected <- matrix(c(0, 5, 5, 5, 5, 0, 0), nrow = 1)

  expect_equal(
    PICsnATAC:::.prepare_pic_counts(dense),
    expected
  )
  expect_equal(
    PICsnATAC:::.prepare_pic_counts(dense, artifact_threshold = Inf),
    matrix(c(0, 5, 5, 5, 5, 5, 5), nrow = 1)
  )

  sparse <- Matrix::Matrix(dense, sparse = TRUE)
  expect_equal(
    as.matrix(PICsnATAC:::.prepare_pic_counts(sparse)),
    expected
  )
  expect_error(
    PICsnATAC:::.prepare_pic_counts(dense, artifact_threshold = 4),
    "at least count_cap"
  )
  expect_error(
    PICsnATAC:::.prepare_pic_counts(matrix(0.5)),
    "whole numbers"
  )
})


test_that("MLE helpers accept dense and sparse matrices consistently", {
  pic_mat <- matrix(c(
    0, 1, 2, 1, 0, 1,
    1, 0, 1, 0, 1, 2
  ), nrow = 2, byrow = TRUE)
  capturing_rates <- rep(0.5, ncol(pic_mat))
  peak_lengths <- c(100, 200)

  dense_loglik <- obs_to_insertion_MLE_obj(
    pic_mat, capturing_rates, peak_lengths
  )
  sparse_loglik <- obs_to_insertion_MLE_obj(
    Matrix::Matrix(pic_mat, sparse = TRUE),
    capturing_rates,
    peak_lengths
  )
  dense_rate <- obs_to_insertion_MLE_lam(
    pic_mat, capturing_rates, peak_lengths
  )
  sparse_rate <- obs_to_insertion_MLE_lam(
    Matrix::Matrix(pic_mat, sparse = TRUE),
    capturing_rates,
    peak_lengths
  )

  expect_equal(sparse_loglik, dense_loglik)
  expect_equal(sparse_rate, dense_rate)
  expect_true(all(is.finite(dense_loglik)))
  expect_true(all(dense_rate >= 0.01 & dense_rate <= 20))
})


test_that("DAR_by_LRT applies the configured artifact threshold", {
  pic_mat <- matrix(c(
    0, 1, 11, 1, 0, 1,
    1, 0, 2, 0, 1, 2
  ), nrow = 2, byrow = TRUE)
  rownames(pic_mat) <- c("chr1:1-101", "chr-2:201-401")
  labels <- rep(c("a", "b"), each = 3)
  capturing_rates <- rep(0.5, ncol(pic_mat))

  default_result <- DAR_by_LRT(
    pic_mat,
    capturing_rates,
    labels,
    plen = c(100, 200)
  )
  manually_prepared <- PICsnATAC:::.prepare_pic_counts(pic_mat)
  manual_result <- DAR_by_LRT(
    manually_prepared,
    capturing_rates,
    labels,
    plen = c(100, 200),
    artifact_threshold = Inf
  )

  expect_equal(default_result, manual_result)
  expect_identical(names(default_result), rownames(pic_mat))
  expect_true(all(default_result >= 0 & default_result <= 1))
  expect_error(
    DAR_by_LRT(pic_mat, capturing_rates, rep("a", ncol(pic_mat)),
      plen = c(100, 200)
    ),
    "exactly two"
  )
})


test_that("moment estimation returns a named matrix", {
  pic_mat <- matrix(c(
    0, 1, 1, 0,
    1, 0, 2, 1
  ), nrow = 2, byrow = TRUE)
  rownames(pic_mat) <- c("chr-foo:1-201", "chr2:100-400")

  result <- obs_to_insertion_ME(
    pic_mat = pic_mat,
    capturing_rates = rep(0.5, 4),
    cell_type_labels = c("a", "a", "b", "b"),
    cap_insertion = 6,
    insertion_rates = c(0.1, 0.2, 0.3),
    peak_lengths = c(200, 300)
  )

  expect_identical(dim(result), c(2L, 2L))
  expect_identical(rownames(result), rownames(pic_mat))
  expect_identical(colnames(result), c("a", "b"))
  expect_true(all(result %in% c(0.1, 0.2, 0.3)))
})
