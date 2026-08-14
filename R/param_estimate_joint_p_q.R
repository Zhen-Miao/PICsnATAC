## param_estimate_joint_p_q


#' Get region (peak) by cell type matrix and per cell capturing rate
#' @description For a snATAC-seq (binary) dataset, compute the peak-specific
#'  open probability
#'  and cell-specific capturing rates
#'
#' @param cell_type_set A vector containing all cell types
#' @param r_by_c Input region (peak) by cell matrix. Nonzero entries are
#'  binarized internally.
#' @param cell_type_labels A vector containing cell type labels
#' @param n_features_per_cell The number of features in the matrix,
#'  defaulting to `nrow(r_by_c)`. If supplied, it must match that value.
#' @param p_acc The accuracy of p, default specified as 0.0005
#' @param q_acc The accuracy of q, default specified as 0.0005
#' @param n_max_iter The maximum iteration, default = 800
#' @param verbose Whether to output information on processing status
#'
#' @return A list with two elements, \itemize{
#'   \item p_by_t Peak by cell type matrix, each element represents
#'    the open probability of the peak in the corresponding cell type
#'   \item q_vec A vector of cell-specific capturing rate
#' }
#' @export
#'
get_r_by_ct_mat_pq <- function(cell_type_set,
                               r_by_c,
                               cell_type_labels,
                               n_features_per_cell = nrow(r_by_c),
                               p_acc = 0.0005,
                               q_acc = 0.0005,
                               n_max_iter = 800,
                               verbose = TRUE) {
  if (!is.matrix(r_by_c) && !methods::is(r_by_c, "Matrix")) {
    stop("r_by_c must be a matrix or Matrix object", call. = FALSE)
  }
  if (nrow(r_by_c) == 0L || ncol(r_by_c) == 0L) {
    stop("r_by_c must have at least one row and one column", call. = FALSE)
  }
  ## require cell names provided
  if (is.null(colnames(r_by_c))) {
    stop("r_by_c must have cell names as column names", call. = FALSE)
  }
  if (anyNA(colnames(r_by_c)) || any(!nzchar(colnames(r_by_c))) ||
      anyDuplicated(colnames(r_by_c))) {
    stop("r_by_c column names must be unique and non-missing", call. = FALSE)
  }
  if (length(cell_type_labels) != ncol(r_by_c) || anyNA(cell_type_labels)) {
    stop("cell_type_labels must contain one non-missing label per cell",
      call. = FALSE
    )
  }
  cell_type_labels <- as.character(cell_type_labels)
  cell_type_set <- as.character(cell_type_set)
  if (length(cell_type_set) == 0L || anyNA(cell_type_set) ||
      any(!nzchar(cell_type_set)) || anyDuplicated(cell_type_set)) {
    stop("cell_type_set must contain unique, non-missing labels", call. = FALSE)
  }
  if (!setequal(cell_type_set, unique(cell_type_labels))) {
    stop("cell_type_set must match the labels present in cell_type_labels",
      call. = FALSE
    )
  }
  if (length(n_features_per_cell) != 1L ||
      n_features_per_cell != nrow(r_by_c)) {
    stop("n_features_per_cell must equal nrow(r_by_c)", call. = FALSE)
  }
  if (length(p_acc) != 1L || !is.finite(p_acc) || p_acc <= 0 ||
      length(q_acc) != 1L || !is.finite(q_acc) || q_acc <= 0) {
    stop("p_acc and q_acc must be positive finite scalars", call. = FALSE)
  }
  if (length(n_max_iter) != 1L || !is.finite(n_max_iter) ||
      n_max_iter < 2L || n_max_iter != as.integer(n_max_iter)) {
    stop("n_max_iter must be an integer of at least 2", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  if (anyNA(r_by_c) || any(!is.finite(r_by_c)) || any(r_by_c < 0)) {
    stop("r_by_c must contain finite, non-negative values", call. = FALSE)
  }

  ## save data to matrix
  itermat_q_by_type <- numeric(ncol(r_by_c))
  names(itermat_q_by_type) <- colnames(r_by_c)
  itermat_p_by_type <- matrix(
    nrow = n_features_per_cell,
    ncol = length(cell_type_set),
    dimnames = list(rownames(r_by_c), cell_type_set)
  )

  ## make the matrix binary
  r_by_c <- Matrix::Matrix(r_by_c, sparse = TRUE)
  r_by_c <- methods::as(r_by_c, "dMatrix")
  r_by_c <- Matrix::drop0(r_by_c)
  r_by_c@x[] <- 1

  ## for each cell type
  for (gg in cell_type_set) {
    r_by_c_sub <- r_by_c[, cell_type_labels == gg, drop = FALSE]
    n_cell_sub <- ncol(r_by_c_sub)
    cell_names_sub <- colnames(r_by_c_sub)
    n_reads_in_cell <- Matrix::colSums(r_by_c_sub)
    n_reads_in_region <- Matrix::rowSums(r_by_c_sub)

    max_reads <- max(n_reads_in_cell)
    if (max_reads == 0) {
      stop(
        "Capturing rates cannot be estimated for cell type '", gg,
        "' because all of its cells have zero accessible peaks",
        call. = FALSE
      )
    }

    ## Starting values without known missing-rate information.
    q_current <- n_reads_in_cell / max_reads
    p_current <- n_reads_in_region / n_cell_sub
    diff1 <- Inf
    diff2 <- Inf
    converged <- FALSE

    for (numiters in seq_len(n_max_iter - 1L)) {
      sum_p <- sum(p_current)
      sum_q <- sum(q_current)
      if (sum_p <= 0 || sum_q <= 0) {
        stop("The p/q iteration encountered a zero denominator", call. = FALSE)
      }

      ## First step -- estimate the missing rate from the open probabilities
      q_new <- pmin(n_reads_in_cell / sum_p, 0.999)

      ## Second step
      p_new <- pmin(n_reads_in_region / sum_q, 0.999)

      ## Use a stable relative change calculation for zero-valued entries.
      diff1 <- mean(
        abs(p_new - p_current) / pmax(abs(p_new), .Machine$double.eps)
      )
      diff2 <- mean(
        abs(q_new - q_current) / pmax(abs(q_new), .Machine$double.eps)
      )
      p_current <- p_new
      q_current <- q_new

      if (diff1 <= p_acc && diff2 <= q_acc) {
        converged <- TRUE
        break
      }
    }

    itermat_p_by_type[, gg] <- p_current
    itermat_q_by_type[cell_names_sub] <- q_current

    if (!converged) {
      warning(
        "p/q estimation for cell type '", gg,
        "' did not converge within ", n_max_iter, " iterations",
        call. = FALSE
      )
    }

    ## print progress
    if (verbose) {
      message(sprintf(
        "%s completed (p change = %.6g; q change = %.6g)",
        gg, diff1, diff2
      ))
    }
  }
  list(p_by_t = itermat_p_by_type, q_vec = itermat_q_by_type)
}
