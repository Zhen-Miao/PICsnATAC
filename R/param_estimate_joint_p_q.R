## param_estimate_joint_p_q


#' Get region (peak) by cell type matrix and per cell capturing rate
#' @description For a snATAC-seq (binary) dataset, compute the peak-specific open probability
#'  and cell-specific capturing rates
#'
#' @param cell_type_set A vector containing all cell types
#' @param r_by_c Input matrix, region (peak) by cell
#' @param cell_type_labels A vector containing cell type labels
#' @param n_features_per_cell The number of features in the matrix, can be calculated by dim(r_by_c)[1]
#' @param p_acc The accuracy of p, default specified as 0.0005
#' @param q_acc The accuracy of q, default specified as 0.0005
#' @param n_max_iter The maximum iteration, default = 500
#'
#' @return A list with two elements, \itemize{
#'   \item p_by_t Peak by cell type matrix, each element represents the open probability of the peak in the corresponding cell type
#'   \item q_vec A vector of cell-specific capturing rate
#' }
#' @export
#'
#' @examples

get_r_by_ct_mat_pq <- function(
  cell_type_set,
  r_by_c,
  cell_type_labels,
  n_features_per_cell,
  p_acc = 0.0005,
  q_acc = 0.0005,
  n_max_iter = 500
){

  ## require r_by_c to have column names
  if(is.null(colnames(r_by_c))){
    stop('the peak by cell matrix has to have column names')
  }
  ## save arguments
  ## note, with this method, there cannot be any unused argument,
  ## or error will happen
  args_s <- list()
  args_s <- append(args_s, mget(names(formals()),sys.frame(sys.nframe())))

  ## get output with internal functions
  outs <- do.call(get_r_by_ct_mat_pq_int, args = args_s)
  itermat_p_by_type <- outs$itermat_p_by_type
  itermat_q_by_type <- outs$itermat_q_by_type
  print(str(itermat_p_by_type))

  ## remove total na values
  for(gg in cell_type_set){
    mat_to_shrink <- itermat_p_by_type[[gg]]
    mat_to_shrink <- mat_to_shrink[,!is.na(colsums(mat_to_shrink))]
    itermat_p_by_type[[gg]] <- mat_to_shrink

    mat_to_shrink2 <- itermat_q_by_type[[gg]]
    mat_to_shrink2 <- mat_to_shrink2[!is.na(rowsums(mat_to_shrink2)),]
    itermat_q_by_type[[gg]] <- mat_to_shrink2
  }

  ## create new matrix to store the p
  p_by_t_new <- matrix(nrow = n_features_per_cell,ncol = length(cell_type_set))
  colnames(p_by_t_new) <- cell_type_set
  for(gg in cell_type_set){
    p_by_t_new_w0 <- itermat_p_by_type[[gg]][,dim(itermat_p_by_type[[gg]])[2]]
    p_by_t_new_w0[p_by_t_new_w0 < 0.0005] = 0.0005
    p_by_t_new[,gg] <- p_by_t_new_w0
  }

  # save q information
  q_vec_new <- vector(length = 0L)
  for(gg in cell_type_set){
    q_vec_op <- itermat_q_by_type[[gg]][dim(itermat_q_by_type[[gg]])[1],]
    clabel <- colnames(r_by_c)[cell_type_labels == gg]
    names(q_vec_op) <- clabel
    # print(sum(is.na(q_vec_op)))
    q_vec_new <- append(q_vec_new, q_vec_op)
  }
  q_vec_new <-  q_vec_new[colnames(r_by_c)]
  return(list(p_by_t = p_by_t_new, q_vec = q_vec_new))
}

#' Get theoretical snATAC-seq distribution under condition 1
#'
#' @importFrom methods is
#' @import Matrix
#' @param insertion_rate The insertion rate (per base pair)
#' @param peak_length The width of peak
#'
#' @return A vector of length 6, representing the probability of observing 0 to >=5 counts
#' @export
#'
#' @examples
#' @keywords internal
#' @noRd
get_r_by_ct_mat_pq_int <- function(
  cell_type_set,
  r_by_c,
  cell_type_labels,
  n_features_per_cell,
  p_acc = 0.0005,
  q_acc = 0.0005,
  n_max_iter = 500
){
  ## save data to matrix
  itermat_q_by_type <- rep(list(),length = length(cell_type_set))
  names(itermat_q_by_type) <- cell_type_set
  itermat_p_by_type <- rep(list(),length = length(cell_type_set))
  names(itermat_p_by_type) <- cell_type_set

  ## make the matrix binary
  if(!is(r_by_c, "sparseMatrix")){
      as(r_by_c, "sparseMatrix")
  }
  r_by_c@x = rep(r_by_c@x, length = length(r_by_c@x))

  ## for each cell type
  for(gg in cell_type_set  ){
    r_by_c_sub  <- r_by_c[,cell_type_labels == gg]
    n_cell_sub <- dim(r_by_c_sub)[2]
    r_by_c_sub <- as.matrix(r_by_c_sub)
    n_reads_in_cell <- colSums(r_by_c_sub)
    n_reads_in_region <- rowSums(r_by_c_sub)

    ## without the true missing rate information
    itermat_q <- matrix(NA,nrow=n_max_iter,ncol=n_cell_sub)
    itermat_p <- matrix(NA,nrow=n_features_per_cell,ncol=n_max_iter)

    itermat_q[1,] <- n_reads_in_cell/max(n_reads_in_cell)  # starting value for q
    itermat_p[,1] <- n_reads_in_region/n_cell_sub    # starting value for p
    diff1 <- 1
    diff2 <- 1
    numiters <- 1

    while ((diff1 > 0.0005 | diff2 > 0.0005) & numiters < n_max_iter){
      q0 <- itermat_q[numiters,]
      p0 <- itermat_p[,numiters]

      ## First step -- estimate the missing rate from the open probabilities
      q0new <- n_reads_in_cell / sum(p0)
      q0new[q0new > 1] <- 0.999 ## make sure it does not exceed 1

      ## Second step
      p0new <- n_reads_in_region / sum(q0)
      p0new[p0new > 1] <- 0.999 ## make sure it does not exceed 1

      ## record the values
      numiters <- numiters + 1
      itermat_p[,numiters] <- p0new
      itermat_q[numiters,] <- q0new

      # print(numiters)
      diff1 <- sum(abs(itermat_p[,numiters]-itermat_p[,numiters-1])/abs(itermat_p[,numiters]), na.rm = T) / n_features_per_cell
      diff2 <- sum(abs(itermat_q[numiters,]-itermat_q[numiters-1,])/abs(itermat_q[numiters,]), na.rm = T) / n_cell_sub
    }
    itermat_p_by_type[[gg]] <- itermat_p
    itermat_q_by_type[[gg]] <- itermat_q
    print(paste(gg,' completed'))
    print(paste('diff1 = ', diff1))
    print(paste('diff2 = ', diff2))
  }
  return(list(itermat_p_by_type = itermat_p_by_type, itermat_q_by_type = itermat_q_by_type))

}
