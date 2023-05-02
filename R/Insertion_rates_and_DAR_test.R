
#' get theoretical snATAC-seq distribution under condition 1
#'
#' @importFrom stats convolve optimize pchisq
#' @param insertion_rate The insertion rate (per 1000 base pairs)
#' @param peak_length The width of peak
#'
#' @return A vector of length 6, representing the probability of observing 0 to >=5 counts
#' @export
#'
#' @keywords internal
#' @noRd
.get_theoretical_c1 <- function(insertion_rate,
                                peak_length = 500) {
  lambda_ <- insertion_rate * peak_length/1000
  p_W_m <- vector(mode = "numeric", length = 6)
  p_W_m[1] <- exp(-0.5 * lambda_) * (2 - exp(-0.5 * lambda_)) ## prob of 0
  p_W_m[2] <- exp(-0.5 * lambda_) * (lambda_ - 2 + 2 * exp(-0.5 * lambda_))
  p_W_m[3] <- exp(-1 * lambda_) * (exp(0.5 * lambda_) * (lambda_^2 - 4 * lambda_ + 8) - 8) / 4
  p_W_m[4] <- exp(-1 * lambda_) * (exp(0.5 * lambda_) * (lambda_^3 - 6 * lambda_^2 + 24 * lambda_ - 48) + 48) / 24
  p_W_m[5] <- exp(-1 * lambda_) * (exp(0.5 * lambda_) * (lambda_^4 - 8 * lambda_^3 + 48 * lambda_^2 - 192 * lambda_ + 384) + -384) / 192
  p_W_m[6] <- 1 - sum(p_W_m) ## prob of observing 5 or higher

  ## bi-allelic scenario
  p_W_m_conv <- convolve(p_W_m, rev(p_W_m), type = "open")
  p_W_m_conv[6] <- sum(p_W_m_conv[6:11])

  ## return a vector of length 6, representing the probability of observing 0 to >=5 counts
  return(p_W_m_conv[1:6])
}


#' Construct an insertion-rate by mean count matrix
#'
#' @param insertion_rates A vector of insertion rates (per 1000 base pairs)
#' @param peak_lengths A vector containing width of peaks
#'
#' @return An insertion-rate by mean count matrix
#' @export
#'
#' @keywords internal
#' @noRd
.insertion_to_c1 <- function(insertion_rates = (1:2000) * 0.01,
                             peak_lengths = 4:20 * 50) {
  ## initialize the matrix
  theo_mean_c1 <- matrix(
    nrow = length(insertion_rates),
    ncol = length(peak_lengths)
  )

  ## iteratively construct the matrix
  for (iii in seq_along(insertion_rates)) {
    for (jjj in seq_along(peak_lengths)) {
      p_W_s_theo <- .get_theoretical_c1(
        insertion_rate = insertion_rates[iii],
        peak_length = peak_lengths[jjj]
      )
      theo_mean_c1[iii, jjj] <- sum(p_W_s_theo * 0:6)
    }
  }
  rownames(theo_mean_c1) <- insertion_rates
  colnames(theo_mean_c1) <- peak_lengths

  return(theo_mean_c1)
}


#' Compute log loss for a given estimated insertion rate
#'
#' @param est_inser Estimated insertion rate (per 1000 base pairs)
#' @param peak_length The width of peak
#' @param cap_fragment Maximum possible fragment (PIC) counts, default = 5
#' @param capturing_rates A vector of capturing rates in each cell
#' @param obs_pic_vec A vector specifying observed fragment (PIC) counts
#'
#' @return Log loss value
#' @export
#'
#' @keywords internal
#' @noRd
.log_loss_frag <- function(est_inser,
                           peak_length,
                           cap_fragment = 5,
                           capturing_rates,
                           obs_pic_vec) {
  ## make sure cap rate has same length with obs_pic
  # if(length(capturing_rates) != length(obs_pic_vec)){
  #   stop('input capturing rate should have same length with observed vector')
  # }
  n_cells <- length(obs_pic_vec)

  p_W_s_theo <- .get_theoretical_c1(
    insertion_rate = est_inser,
    peak_length = peak_length
  )

  ## consider missing
  lg_p_W_o_t <- vector(mode = "numeric", length = 9)
  all_capturing_rates <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)

  for (c_rate in 1:9) {
    ct <- all_capturing_rates[c_rate]
    pt <- matrix(0, nrow = cap_fragment + 1, ncol = cap_fragment + 1)
    rownames(pt) <- colnames(pt) <- 0:cap_fragment
    for (ot in 0:cap_fragment) {
      for (k in ot:cap_fragment) {
        pt[as.character(ot), as.character(k)] <- choose(k, ot) * ct^ot * (1 - ct)^(k - ot) * p_W_s_theo[k + 1]
      }
    }
    ptv <- rowSums(pt)
    obs_pic_vec_sel <- as.character(obs_pic_vec[capturing_rates == ct])
    probs <- log(ptv[obs_pic_vec_sel])

    ## remove the effect of very large value
    probs[probs < -15] <- -15
    lg_p_W_o_t[c_rate] <- sum(probs)
  }

  ll <- sum(lg_p_W_o_t)
  return(ll)
}


#' Calculate the optimal loss from observed PIC counts
#'
#' @import parallel
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param plen A vector of peak widths
#' @param n_cores A numerical value to specify the number of cores in parallel
#'
#' @return The optimized loss over insertion rates from 0.01 to 20
#' @export
#'
obs_to_insertion_MLE_obj <- function(pic_mat,
                                 capturing_rates,
                                 plen,
                                 n_cores) {
  n_para <- length(plen)
  optim_results <- vector(length = n_para)

  ## setup mltiple cores
  cl <- makeCluster(n_cores)
  clusterExport(cl, c(".log_loss_frag"))

  ## iterations
  optim_results <- mclapply(1:n_para, function(pp) {
    optimize(
      f = .log_loss_frag,
      interval = c(0.01 , 20),
      peak_length = plen[pp],
      capturing_rates = capturing_rates,
      obs_pic_vec = pic_mat[pp, ],
      maximum = T
    )$objective
  }, mc.cores = n_cores)
  optim_results = unlist(optim_results)
  stopCluster(cl)
  return(optim_results)
}


#' Calculate the insertion rate from observed PIC counts
#'
#' @import parallel
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param plen A vector of peak widths
#' @param n_cores A numerical value to specify the number of cores in parallel
#'
#' @return The optimized insertion rate (per 1000 base pairs) over insertion rates
#'      from 0.01 to 20
#' @export
#'
obs_to_insertion_MLE_lam <- function(pic_mat,
                                     capturing_rates,
                                     plen,
                                     n_cores) {
  n_para <- length(plen)
  optim_results <- vector(length = n_para)

  ## setup mltiple cores
  cl <- makeCluster(n_cores)
  clusterExport(cl, c(".log_loss_frag"))

  ## iterations
  optim_results <- mclapply(1:n_para, function(pp) {
    optimize(
      f = .log_loss_frag,
      interval = c(0.01, 20),
      peak_length = plen[pp],
      capturing_rates = capturing_rates,
      obs_pic_vec = pic_mat[pp, ],
      maximum = T
    )$maximum
  }, mc.cores = n_cores)
  optim_results = unlist(optim_results)
  stopCluster(cl)
  return(optim_results)
}

#' Compute p value for DAR test between two cell types
#'
#' @import parallel
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param cell_type_labels A vector specifying cell type labels
#' @param estimation_approach The approach for parameter estimation
#' @param n_cores A numerical value to specify the number of cores in parallel
#' @param plen A vector of peak length
#'
#' @return Log loss value
#' @export
#'
DAR_by_LRT <- function(pic_mat,
                       capturing_rates,
                       cell_type_labels,
                       n_cores,
                       plen = NULL,
                       estimation_approach = "MLE") {
  ## load library
  requireNamespace('parallel')

  ## save some values
  n_pks <- dim(pic_mat)[1]
  ct_uniq <- unique(cell_type_labels)

  ## pk length
  if(is.null(plen)){
    plen <- rep(500, times = n_pks)
  }


  ## cap the counts at 5
  cts <- pic_mat@x
  cts[cts >= 10] <- 0 ## counts greater than 10 are likely artifact
  cts[cts > 5] <- 5
  pic_mat@x <- cts
  rm(cts)

  ## split the matrix by cell types
  pic_mat <- as.matrix(pic_mat)
  pic_mat_1 <- pic_mat[, cell_type_labels == ct_uniq[1]]
  pic_mat_2 <- pic_mat[, cell_type_labels == ct_uniq[2]]

  ## correct for the capturing rate for fast computation
  capturing_rates <- ceiling(capturing_rates * 10) / 10
  capturing_rates[capturing_rates < 0.2] <- 0.2
  capturing_rates[capturing_rates > 0.9] <- 0.99

  capturing_rates_1 <- capturing_rates[cell_type_labels == ct_uniq[1]]
  capturing_rates_2 <- capturing_rates[cell_type_labels == ct_uniq[2]]

  ## MLE
  if (estimation_approach == "MLE") {
    ## likelihood under the null model
    ll_all_mle <- obs_to_insertion_MLE_obj(
      pic_mat = pic_mat,
      capturing_rates = capturing_rates,
      plen = plen,
      n_cores = n_cores
    )
    ## likelihood under the full model (alternative)
    ll_full_1_mle <- obs_to_insertion_MLE_obj(
      pic_mat = pic_mat_1,
      capturing_rates = capturing_rates[cell_type_labels == ct_uniq[1]],
      plen = plen,
      n_cores = n_cores
    )

    ll_full_2_mle <- obs_to_insertion_MLE_obj(
      pic_mat = pic_mat_2,
      capturing_rates = capturing_rates[cell_type_labels == ct_uniq[2]],
      plen = plen,
      n_cores = n_cores
    )

    ## p value is obtained by chi-squared statistics
    p_val <- pchisq(2 * (ll_full_1_mle + ll_full_2_mle - ll_all_mle), df = 1, lower.tail = F)
  }

  return(p_val)
}
