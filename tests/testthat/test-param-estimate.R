test_that("p/q estimation binarizes dense and sparse count matrices", {
  counts <- matrix(c(
    0, 2, 0, 1,
    3, 0, 1, 0,
    1, 1, 0, 0
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("cell-", seq_len(ncol(counts)))
  rownames(counts) <- paste0("peak-", seq_len(nrow(counts)))
  labels <- c("a", "a", "b", "b")

  dense_result <- get_r_by_ct_mat_pq(
    c("a", "b"),
    counts,
    labels,
    verbose = FALSE
  )
  binary_result <- get_r_by_ct_mat_pq(
    c("a", "b"),
    (counts != 0) * 1,
    labels,
    verbose = FALSE
  )
  sparse_result <- get_r_by_ct_mat_pq(
    c("a", "b"),
    Matrix::Matrix(counts, sparse = TRUE),
    labels,
    verbose = FALSE
  )

  expect_equal(dense_result, binary_result)
  expect_equal(sparse_result, binary_result)
  expect_identical(rownames(dense_result$p_by_t), rownames(counts))
  expect_true(all(is.finite(dense_result$p_by_t)))
  expect_true(all(dense_result$p_by_t >= 0 & dense_result$p_by_t <= 1))
  expect_true(all(dense_result$q_vec >= 0 & dense_result$q_vec <= 1))
})


test_that("p/q estimation rejects inconsistent inputs", {
  counts <- matrix(c(1, 0, 0, 1), nrow = 2)
  colnames(counts) <- c("cell-a", "cell-b")

  expect_error(
    get_r_by_ct_mat_pq("a", counts, c("a", "b"), verbose = FALSE),
    "must match"
  )
  expect_error(
    get_r_by_ct_mat_pq(c("a", "b"), counts, c("a"), verbose = FALSE),
    "one non-missing label per cell"
  )
  expect_error(
    get_r_by_ct_mat_pq(
      c("a", "b"), counts, c("a", "b"),
      n_features_per_cell = 3,
      verbose = FALSE
    ),
    "must equal"
  )
})
