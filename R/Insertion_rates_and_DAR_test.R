
#' Title get theoretical snATAC-seq distribution under condition 1
#'
#' @param insertion_rate The insertion rate (per base pair)
#' @param peak_length The width of peak
#'
#' @return A vector of length 6, representing the probability of observing 0 to >=5 counts
#' @export
#'
#' @examples
#' @keywords internal
.get_theoretical_c1 <- function(insertion_rate,
                                peak_length = 500){
  lambda_ = insertion_rate * peak_length
  p_W_m = vector(mode = 'numeric',length = 6)
  p_W_m[1] = exp(-0.5*lambda_) * (2-exp(-0.5*lambda_)) ## prob of 0
  p_W_m[2] = exp(-0.5*lambda_) * (lambda_ - 2 + 2 * exp(-0.5*lambda_))
  p_W_m[3] = exp(-1*lambda_)*( exp(0.5*lambda_) *(lambda_^2-4*lambda_+8) -8  )/4
  p_W_m[4] = exp(-1*lambda_)*( exp(0.5*lambda_) *(lambda_^3 -6 * lambda_^2 +24*lambda_-48) +48  )/24
  p_W_m[5] = exp(-1*lambda_)*( exp(0.5*lambda_) *(lambda_^4 -8 * lambda_^3 +48*lambda_^2 -192 * lambda_+384) +-384  )/192
  p_W_m[6] = 1 - sum(p_W_m) ## prob of observing 5 or higher

  ## bi-allelic scenario
  p_W_m_conv <- convolve(p_W_m, rev(p_W_m), type = 'open')
  p_W_m_conv[6] = sum(p_W_m_conv[6:11])

  ## return a vector of length 6, representing the probability of observing 0 to >=5 counts
  return(p_W_m_conv[1:6])
}


#' Title Construct an insertion-rate by mean count matrix
#'
#' @param insertion_rates A vector of insertion rates (per base pair)
#' @param peak_lengths A vector containing width of peaks
#'
#' @return An insertion-rate by mean count matrix
#' @export
#'
#' @examples
#' @keywords internal
.insertion_to_c1 <-  function(
    insertion_rates = (1:2000)*0.01*0.001,
    peak_lengths = 4:20*50
){
  ## initialize the matrix
  theo_mean_c1 <-  matrix(nrow = length(insertion_rates),
                          ncol = length(peak_lengths))

  ## iteratively construct the matrix
  for(iii in seq_along(insertion_rates)){
    for(jjj in seq_along(peak_lengths)){
      p_W_s_theo <- .get_theoretical_c1(
        insertion_rate = insertion_rates[iii],
        peak_length = peak_lengths[jjj])
      theo_mean_c1[iii,jjj] = sum(p_W_s_theo * 0:6)
    }
  }
  rownames(theo_mean_c1) <- insertion_rates
  colnames(theo_mean_c1) <- peak_lengths

  return(theo_mean_c1)
}


#' Title Compute log loss for a given estimated insertion rate
#'
#' @param est_inser Estimated insertion rate
#' @param peak_length The width of peak
#' @param cap_fragment Maximum possible fragment (PIC) counts, default = 5
#' @param capturing_rates A vector of capturing rates in each cell
#' @param obs_pic_vec A vector specifying observed fragment (PIC) counts
#'
#' @return Log loss value
#' @export
#'
#' @examples
#' @keywords internal
.log_loss_frag <- function(
    est_inser,
    peak_length,
    cap_fragment = 5,
    capturing_rates,
    obs_pic_vec
){
  ## make sure cap rate has same length with obs_pic
  # if(length(capturing_rates) != length(obs_pic_vec)){
  #   stop('input capturing rate should have same length with observed vector')
  # }
  n_cells <- length(obs_pic_vec)

  p_W_s_theo <- .get_theoretical_c1(insertion_rate = est_inser,
                                    peak_length = peak_length)

  ## consider missing
  lg_p_W_o_t <- vector(mode = 'numeric',length = 9)
  all_capturing_rates = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)

  for(c_rate in 1:9){
    ct = all_capturing_rates[c_rate]
    pt = matrix(0, nrow = cap_fragment+1, ncol = cap_fragment+1)
    rownames(pt) <- colnames(pt) <- 0:cap_fragment
    for(ot in 0:cap_fragment){
      for(k in ot:cap_fragment ){
        pt[as.character(ot),as.character(k)] <- choose(k,ot) * ct^ot * (1-ct)^(k-ot) * p_W_s_theo[k+1]
      }
    }
    ptv <- rowSums(pt)
    obs_pic_vec_sel = as.character(obs_pic_vec[capturing_rates == ct])
    probs <- log(ptv[obs_pic_vec_sel])

    ## remove the effect of very large value
    probs[probs < -15] <- -15
    lg_p_W_o_t[c_rate] = sum(probs)
  }

  ll = sum(lg_p_W_o_t)
  return(ll)
}


#' Title Calculate the insertion rate from observed PIC counts
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param plen A vector of peak widths
#'
#' @return The optimized insertion rate
#' @export
#'
#' @examples
obs_to_insertion_MLE <- function(
    pic_mat,
    capturing_rates,
    plen
){
  n_para = length(plen)
  optim_results = vector(length = n_para)

  ## to-do: use mclapply to speed up the process
  for(pp in 1:n_para){
    optim_results[pp] = optimize(f = .log_loss_frag,
                                 interval = c(0.01*0.001, 0.02),
                                 peak_length = plen[pp],
                                 capturing_rates = capturing_rates,
                                 obs_pic_vec = pic_mat[pp,],
                                 maximum = T)$objective
  }
  return(optim_results)
}

#' Title Compute p value for DAR test between two cell types
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param cell_type_labels A vector specifying cell type labels
#'
#' @return Log loss value
#' @export
#'
#' @examples
DAR_by_LRT <- function(
    pic_mat,
    capturing_rates,
    cell_type_labels
){
  ## save some values
  n_pks <- dim(pic_mat)[1]
  ct_uniq <- unique(cell_type_labels)

  ## pk length
  plen = rep(500, times = n_pks)

  pic_mat <- as.matrix(pic_mat)
  pic_mat_1 <- pic_mat[,cell_type_labels == ct_uniq[1]]
  pic_mat_2 <- pic_mat[,cell_type_labels == ct_uniq[2]]

  ## correct for the capturing rate for fast computation
  capturing_rates = ceiling(capturing_rates *10) / 10
  capturing_rates[capturing_rates < 0.2] <- 0.2
  capturing_rates[capturing_rates > 0.9] <- 0.99

  capturing_rates_1 <- capturing_rates[cell_type_labels == ct_uniq[1]]
  capturing_rates_2 <- capturing_rates[cell_type_labels == ct_uniq[2]]

  ll_all_mle <- obs_to_insertion_MLE(pic_mat = pic_mat,
                                       capturing_rates = capturing_rates,
                                       plen = plen)

  ll_full_1_mle <- obs_to_insertion_MLE(pic_mat = pic_mat_1,
                                          capturing_rates = capturing_rates[cell_type_labels == ct_uniq[1]],
                                          plen = plen)

  ll_full_2_mle <- obs_to_insertion_MLE(pic_mat = pic_mat_2,
                                          capturing_rates = capturing_rates[cell_type_labels == ct_uniq[2]],
                                          plen = plen)

  p_val <- pchisq(2*(ll_full_1_mle+ll_full_2_mle - ll_all_mle), df = 1,lower.tail = F)
  return(p_val)
}



