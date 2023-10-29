## param_estimate_joint_p_q


#' Get region (peak) by cell type matrix and per cell capturing rate
#' @description For a snATAC-seq (binary) dataset, compute the peak-specific
#'  open probability
#'  and cell-specific capturing rates
#'
#' @importFrom methods is as
#' @param cell_type_set A vector containing all cell types
#' @param r_by_c Input matrix, region (peak) by cell
#' @param cell_type_labels A vector containing cell type labels
#' @param n_features_per_cell The number of features in the matrix,
#'  can be calculated by nrow(r_by_c)
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
                               n_features_per_cell,
                               p_acc = 0.0005,
                               q_acc = 0.0005,
                               n_max_iter = 800,
                               verbose = TRUE) {
  ## require cell names provided
  if (is.null(colnames(r_by_c))) {
    stop("the peak by cell matrix has to have column names")
  }

  ## save data to matrix
  itermat_q_by_type <- vector(length = dim(r_by_c)[2])
  names(itermat_q_by_type) <- colnames(r_by_c)
  itermat_p_by_type <- matrix(
    nrow = n_features_per_cell,
    ncol = length(cell_type_set)
  )
  colnames(itermat_p_by_type) <- cell_type_set

  p_by_t_new <- matrix(nrow = n_features_per_cell, ncol = length(cell_type_set))
  colnames(p_by_t_new) <- cell_type_set


  ## make the matrix binary
  if (!is(r_by_c, "sparseMatrix")) {
    as(r_by_c, "sparseMatrix")
  }
  r_by_c@x <- rep(r_by_c@x, length = length(r_by_c@x))

  ## for each cell type
  for (gg in cell_type_set) {
    r_by_c_sub <- r_by_c[, cell_type_labels == gg]
    n_cell_sub <- dim(r_by_c_sub)[2]
    cell_names_sub <- colnames(r_by_c_sub)
    n_reads_in_cell <- colSums(r_by_c_sub)
    n_reads_in_region <- rowSums(r_by_c_sub)

    ## without the true missing rate information
    itermat_q <- matrix(NA, nrow = n_max_iter, ncol = n_cell_sub)
    itermat_p <- matrix(NA, nrow = n_features_per_cell, ncol = n_max_iter)

    itermat_q[1, ] <- n_reads_in_cell / max(n_reads_in_cell) # starting for q
    itermat_p[, 1] <- n_reads_in_region / n_cell_sub # starting value for p
    diff1 <- 1
    diff2 <- 1
    numiters <- 1

    while ((diff1 > p_acc || diff2 > q_acc) && numiters < n_max_iter) {
      q0 <- itermat_q[numiters, ]
      p0 <- itermat_p[, numiters]

      ## First step -- estimate the missing rate from the open probabilities
      q0new <- n_reads_in_cell / sum(p0)
      q0new[q0new > 1] <- 0.999 ## make sure it does not exceed 1

      ## Second step
      p0new <- n_reads_in_region / sum(q0)
      p0new[p0new > 1] <- 0.999 ## make sure it does not exceed 1

      ## record the values
      numiters <- numiters + 1
      itermat_p[, numiters] <- p0new
      itermat_q[numiters, ] <- q0new

      diff1 <- sum(abs(itermat_p[, numiters] - itermat_p[, numiters - 1]) /
                     abs(itermat_p[, numiters]), na.rm = TRUE) /
        n_features_per_cell
      diff2 <- sum(abs(itermat_q[numiters, ] - itermat_q[numiters - 1, ]) /
                     abs(itermat_q[numiters, ]), na.rm = TRUE) / n_cell_sub
    }

    ## remove columns that contain na values
    mat_to_shrink <- itermat_p
    mat_to_shrink <- mat_to_shrink[, !is.na(colSums(mat_to_shrink))]
    itermat_p_by_type[, gg] <- mat_to_shrink[, dim(mat_to_shrink)[2]]

    mat_to_shrink2 <- itermat_q
    mat_to_shrink2 <- mat_to_shrink2[!is.na(rowSums(mat_to_shrink2)), ]
    itermat_q_by_type[cell_names_sub] <- mat_to_shrink2[dim(mat_to_shrink2)[1],]

    ## print progress
    if (verbose) {
      print(paste(gg, " completed"))
      print(paste("diff1 = ", diff1))
      print(paste("diff2 = ", diff2))
    }
  }
  return(list(p_by_t = itermat_p_by_type, q_vec = itermat_q_by_type))
}
