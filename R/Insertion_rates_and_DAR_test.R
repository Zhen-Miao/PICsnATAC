
#' Title get theoretical snATAC-seq distribution under condition 1
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

#' Title get_theoretical_c12
#'
#' @param insertion_rate A single value of insertion rate
#' @param peak_length A single vale of peak length, default 1000
#' @param pad_length Length of pad region, i.e., flanking regions of the peak
#'  where insertion can also happen
#' @param dirs Directions of insertion, for standard ATAC experiment (default),
#'  it is c(-1,1)
#' @param count_type output count type, default is 'PIC'
#' @param min_frag_length minimum fragment length, default 25
#' @param max_frag_length maximum fragment length, default 600
#' @param cap_insertion number of insertions capped at this value, default 20
#'
#' return A vector of probability corresponding to different number of fragments
#' #' @noRd
#'
.get_theoretical_c12 <- function(
    insertion_rate,
    peak_length = 1000,
    pad_length = 0,
    dirs = c(-1,1),
    count_type = 'PIC',
    min_frag_length = 25,
    max_frag_length = 600,
    cap_insertion = 20 ## assume there will not be more than 20 insertions in
    ## one peak region
){
  lambda_ = insertion_rate * peak_length / 1000
  cap_fragment = floor(cap_insertion/2)
  p_W_m = vector(length = cap_fragment)

  for(mi in seq_along(p_W_m)){
    all_index_n = (mi+1):cap_insertion
    p_W_m_n = vector(length = length(all_index_n))
    for(ni in seq_along(all_index_n)){
      index_n = all_index_n[ni]
      p_W_m_n_ni = ( choose(index_n-1, mi) * 0.5^(index_n-1) * lambda_^index_n *
                       exp(-1 * lambda_) ) / (factorial(index_n))
      if(is.na(p_W_m_n_ni)){
        print(insertion_rate)
        # browser()
      }else if(length(p_W_m_n_ni) != 1){
        print(insertion_rate)
        # browser()
      }
      p_W_m_n[ni] = p_W_m_n_ni
    }
    p_W_m[mi] =  sum(p_W_m_n)
  }
  # print(sum(p_W_m * (1:length(p_W_m))))

  p_W_s = vector(length = cap_fragment)

  exp_capture = (exp(-1 * min_frag_length * insertion_rate / 1000) -
                   exp(-1 * max_frag_length * insertion_rate / 1000)) /
    (1 - exp(-1*peak_length * insertion_rate / 1000))
  # print(exp_capture)

  for(ki in seq_along(p_W_s)){
    all_index_m = ki:cap_fragment
    p_W_s_k = vector(length = length(all_index_m))
    for(mii in seq_along(p_W_s_k)){
      index_m = all_index_m[mii]
      p_W_s_k[mii] = choose(index_m, ki) * exp_capture^ki *
        (1- exp_capture)^(index_m - ki) * p_W_m[index_m]
    }
    p_W_s[ki] = sum(p_W_s_k)
  }

  ## bi-allelic scenario, we convolve it with itself
  p_W_s <-  c(1 - sum(p_W_s), p_W_s)
  p_W_s_conv <- convolve(p_W_s, rev(p_W_s), type = 'open')
  return(p_W_s_conv[2:(cap_fragment+1)])
}



#' given insertion rate, calculate theoretical distribution under condition 1
#'  and 2
#'
#' @param insertion_rates A vector of insertion rates,
#'  note that the insertion rate is per 1000 base, so the default range
#'  is from 0.001 to 20
#' @param peak_lengths A vector of peak length ranges, default is
#'  from 200 to 1000
#' @param pad_length how many bp outside of peaks can insertions happen,
#'  default = 0
#' @param dirs directions. default (-1, 1)
#' @param count_type what type of count should we estimate, default 'fragment'
#' @param min_frag_length minimum of fragment length, default 25
#' @param max_frag_length maximum of fragment length, default 600
#' @param cap_insertion maximum number of insertions that can simultaneously
#'  happen within one peak region
#'
#' @return A insertion rate by peak length matrix, and each
#' element of the matrix corresponds to the number of expected fragments
#' #' @noRd
.insertion_to_c12 <-  function(
    insertion_rates = (1:2000)*0.01,
    peak_lengths = 4:20*50,
    pad_length = 0,
    dirs = c(-1,1),
    count_type = 'fragment',
    min_frag_length = 25,
    max_frag_length = 600,
    cap_insertion = 20
){
  ## initialize the matrix
  theo_mean_c12 <-  matrix(nrow = length(insertion_rates),
                           ncol = length(peak_lengths))

  for(iii in seq_along(insertion_rates)){
    for(jjj in seq_along(peak_lengths)){
      p_W_s_theo <- .get_theoretical_c12(
        insertion_rate = insertion_rates[iii],
        peak_length = peak_lengths[jjj],
        pad_length = pad_length,
        dirs = dirs,
        count_type = count_type,
        min_frag_length = min_frag_length,
        max_frag_length = max_frag_length,
        cap_insertion = cap_insertion)
      theo_mean_c12[iii,jjj] = sum(p_W_s_theo * 1:length(p_W_s_theo))
    }
  }
  rownames(theo_mean_c12) <- insertion_rates
  colnames(theo_mean_c12) <- peak_lengths

  return(theo_mean_c12)
}



#' Title Compute log loss for a given estimated insertion rate
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
.log_loss_frag_c1 <- function(est_inser,
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


#' Title Calculate the optimal loss from observed PIC counts
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

  # ## setup mltiple cores
  # cl <- makeCluster(n_cores)
  # clusterExport(cl, c(".log_loss_frag_c1"))

  ## iterations
  optim_results <- mclapply(1:n_para, function(pp) {
    optimize(
      f = .log_loss_frag_c1,
      interval = c(0.01 , 20),
      peak_length = plen[pp],
      capturing_rates = capturing_rates,
      obs_pic_vec = pic_mat[pp, ],
      maximum = T
    )$objective
  }, mc.cores = n_cores)

  optim_results = unlist(optim_results)

  # stopCluster(cl)
  return(optim_results)
}


#' Title Moment estimator for insertion rates from observed values
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param cell_type_labels A vector of cell type lables for each cell
#' @param cap_insertion The maximum number of insertions in a peak region.
#' @param insertion_rates The range of insertion rates (per 1000 bp) to be considered.
#' @param peak_lengths A vector of peak lengths
#' @param min_frag_length The value for the s1 hyperparameter in the ssPoisson
#'   distribution, this stands for the minimum fragment length requirement such
#'   that the fragment can be amplifiable and mappable to genome. Default = 25
#' @param max_frag_length The value for the s2 hyperparameter in the ssPoisson
#'   distribution, this stands for the max fragment length requirement such
#'   that the fragment can be amplifiable. Default = 600

#'
#' @return A matrix of estimated insertion rates
#' @export
#'
obs_to_insertion_ME <- function(
    pic_mat,
    capturing_rates,
    cell_type_labels,
    cap_insertion = 20,
    min_frag_length = 25,
    max_frag_length = 600,
    insertion_rates = (1:2000)*0.01 ,
    peak_lengths = 4:20*50
){
  ## get the matrix
  c12_mat <- .insertion_to_c12(
    insertion_rates = insertion_rates,
    peak_lengths = peak_lengths,
    min_frag_length = min_frag_length,
    max_frag_length = max_frag_length
  )

  ## ------- group peaks into different
  ## get peak length
  pks <- data.frame(ranges = rownames(pic_mat))
  pks <- tidyr::separate(pks, col = "ranges", sep = ':|-',
                         into = c("chr", "start", "end"))
  plen <- as.numeric(pks$end) - as.numeric(pks$start )

  ## assign group
  plen_gap <- peak_lengths[2] - peak_lengths[1]
  plen_group <- ceiling((plen - peak_lengths[1]) / plen_gap) +1
  plen_group[plen_group > length(peak_lengths)] <- length(peak_lengths)
  # names(plen_group) <- rownames(pic_mat)

  ## to prevent the low acc cells, adjust capturing rate
  capturing_rates[capturing_rates < 0.2] <- 0.2

  ## do this for each cell type
  cell_types <- unique(cell_type_labels)
  plen_group <- as.character(plen_group)
  pgroups <- unique(plen_group)

  est_inser_all <- rep(list(), length = length(cell_types))
  names(est_inser_all) <- cell_types

  for(ct in cell_types){
    sel_cell <- cell_type_labels == ct
    pic_sub <- pic_mat[,sel_cell]
    w_bar <- rowSums(pic_sub %*% diag(1 / capturing_rates[sel_cell]))
    w_bar <- w_bar / sum(sel_cell)
    names(w_bar) <- rownames(pic_sub)

    # save output in a list
    est_inser <- rep(list(), length = length(pgroups))
    names(est_inser) <- pgroups
    for(pg in pgroups){
      c12_sub <- c12_mat[,as.numeric(pg)]
      w_bar_sub <- w_bar[plen_group == pg]
      lkup <- outer(w_bar_sub,c12_sub, '-' )
      lkup <- abs(lkup)

      # for each row, look up the smallest value
      est_inser_ind <- apply(lkup, 1, which.min)
      est_inser[[pg]] <- insertion_rates[est_inser_ind]
      names(est_inser[[pg]]) <- names(w_bar_sub)
    }
    est_inser <- est_inser[lengths(est_inser) != 0]
    est_inser <- unname(est_inser)
    est_inser_all[[ct]] <- do.call('c',args = est_inser)
    est_inser_all[[ct]] <- est_inser_all[[ct]][rownames(pic_mat)]

  }

  est_inser_all_matrix <- do.call(cbind,est_inser_all)

  return(est_inser_all_matrix)
}


#' Title
#'
#' @param est_inser Estimated insertion rate
#' @param peak_length A numeric value of peak length
#' @param min_frag_length The value for the s1 hyperparameter in the ssPoisson
#'   distribution, this stands for the minimum fragment length requirement such
#'   that the fragment can be amplifiable and mappable to genome. Default = 25
#' @param max_frag_length The value for the s2 hyperparameter in the ssPoisson
#'   distribution, this stands for the max fragment length requirement such
#'   that the fragment can be amplifiable. Default = 600
#' @param cap_fragment The maximum number of fragments. Values beyond this will
#'   be treated as outlier and set to this value. Default = 10
#' @param capturing_rates The value of capturing probability for each cell
#' @param obs_pic_vec Observed PIC count vector for the peak
#'
#' @return Log loss value
#' @export
#'
#' @keywords internal
#' @noRd
.log_loss_frag_ME <- function(
    est_inser,
    peak_length,
    min_frag_length = 25,
    max_frag_length = 600,
    cap_fragment = 10,
    capturing_rates,
    obs_pic_vec
){
  ## make sure cap rate has same length with obs_pic
  # if(length(capturing_rates) != length(obs_pic_vec)){
  #   stop('input capturing rate should have same length with observed vector')
  # }
  n_cells <- length(obs_pic_vec)



  p_W_s_theo <- .get_theoretical_c12(insertion_rate = est_inser,
                                     peak_length = peak_length,
                                     min_frag_length = min_frag_length,
                                     max_frag_length = max_frag_length)

  ## include p(W_s = 0)
  p_W_s_0 = 1 - sum(p_W_s_theo)
  # if(p_W_s_0 < 0 ){
  #   stop('negative possibility of getting count 0, check the model setting')
  # }
  p_W_s_theo <- c(p_W_s_0, p_W_s_theo)
  lg_p_W_o_t <- vector(mode = 'numeric',length = n_cells)
  for(cell_i in 1:n_cells){
    ot <- obs_pic_vec[cell_i]
    ct <- capturing_rates[cell_i]
    pt <- vector(mode = 'numeric',length = cap_fragment-ot+1)
    names(pt) <- ot:cap_fragment
    for(k in ot:cap_fragment ){
      pt[as.character(k)] <- choose(k,ot) * ct^ot * (1-ct)^(k-ot) * p_W_s_theo[k+1]
    }
    lg_p_W_o_t[cell_i] = log(sum(pt))
  }
  ## remove the effect of very large value
  lg_p_W_o_t[lg_p_W_o_t < -15] <- -15
  ll = sum(lg_p_W_o_t)
  return(ll)
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

  # ## setup mltiple cores
  # cl <- makeCluster(n_cores)
  # clusterExport(cl, c(".log_loss_frag_c1"))

  ## iterations
  optim_results <- mclapply(1:n_para, function(pp) {
    optimize(
      f = .log_loss_frag_c1,
      interval = c(0.01, 20),
      peak_length = plen[pp],
      capturing_rates = capturing_rates,
      obs_pic_vec = pic_mat[pp, ],
      maximum = T
    )$maximum
  }, mc.cores = n_cores)
  optim_results = unlist(optim_results)
  # stopCluster(cl)
  return(optim_results)
}

#' Compute p value for DAR test between two cell types
#'
#' @import parallel
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param cell_type_labels A vector specifying cell type labels
#' @param estimation_approach The approach for parameter estimation, either 'MLE'
#'   for condition 1 or 'ME' for condition 1+2. The 'MLE' approach is more accurate
#'   and usually it has a higher power, but it ignores the size filtering step
#'   in snATAC-seq data generation. Default is 'MLE'
#' @param n_cores A numerical value to specify the number of cores in parallel,
#'   default = 1
#' @param plen A vector of peak length
#' @param min_frag_length The value for the s1 hyperparameter in the ssPoisson
#'   distribution, this stands for the minimum fragment length requirement such
#'   that the fragment can be amplifiable and mappable to genome. Default = 25
#' @param max_frag_length The value for the s2 hyperparameter in the ssPoisson
#'   distribution, this stands for the max fragment length requirement such
#'   that the fragment can be amplifiable. Default = 600
#'
#' @return Log loss value
#' @export
#'
DAR_by_LRT <- function(pic_mat,
                       capturing_rates,
                       cell_type_labels,
                       n_cores = 1,
                       plen = NULL,
                       min_frag_length = 25,
                       max_frag_length = 600,
                       estimation_approach = "MLE") {
  ## load library
  requireNamespace('parallel')

  ## check input
  if(length(capturing_rates) != length(cell_type_labels)|
     dim(pic_mat)[2] != length(capturing_rates)){
    stop('Number of cells do not match among inputs')
  }

  if(!(estimation_approach %in% c('MLE', 'ME'))){
    stop('Please choose estimation_approach from MLE or ME')
  }

  if(nrow(pic_mat) * ncol(pic_mat) >= (2^31-1)){
    print('pic_mat too large, please slice the matrix by peaks and run the test')
    print('i.e., each time, test for different set of peaks ')
    stop('pic_mat too large')
  }

  ## save some values
  n_pks <- dim(pic_mat)[1]
  ct_uniq <- unique(cell_type_labels)

  ## pk length
  if(is.null(plen)){
    pks <- data.frame(ranges = rownames(pic_mat))
    pks <- tidyr::separate(pks, col = "ranges", sep = ':|-',
                           into = c("chr", "start", "end"))
    plen <- as.numeric(pks$end) - as.numeric(pks$start)
    ## we only need to group peak length by n*100 bp
    plen = ceiling(plen/100)*100
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
  }else if(estimation_approach == "ME"){

    ## calculate estimated lambda for null hypothesis
    lamb_all <- obs_to_insertion_ME(
      pic_mat = pic_mat,
      capturing_rates = capturing_rates,
      min_frag_length = min_frag_length,
      max_frag_length = max_frag_length,
      cell_type_labels = rep('A', length = length(capturing_rates))
    )
    lamb_all <- lamb_all[,1]

    if(anyNA(lamb_all)){
      print('NA detected in lamb_all')
      print(which(is.na(lamb_all)))
    }

    ## calculate estimated lambda for alternative hypothesis
    lamb_full <- obs_to_insertion_ME(
      pic_mat = pic_mat,
      capturing_rates = capturing_rates,
      min_frag_length = min_frag_length,
      max_frag_length = max_frag_length,
      cell_type_labels = cell_type_labels
    )
    lamb_full_1 <- lamb_full[,1]
    lamb_full_2 <- lamb_full[,2]


    if(anyNA(lamb_full_1)){
      print('NA detected in lamb_full_1')
      print(which(is.na(lamb_full_1)))
    }else if(anyNA(lamb_full_2)){
      print('NA detected in lamb_full_2')
      print(which(is.na(lamb_full_2)))
    }


    capturing_rates_1 <- capturing_rates[cell_type_labels == ct_uniq[1]]
    capturing_rates_2 <- capturing_rates[cell_type_labels == ct_uniq[2]]
    pic_mat_1 <- pic_mat[,cell_type_labels == ct_uniq[1]]
    pic_mat_2 <- pic_mat[,cell_type_labels == ct_uniq[2]]

    ## test for each peak
    ll_null <- vector(mode = 'numeric', length = n_pks)
    p_val <- ll_full <- ll_null


    for(pki in 1:n_pks){
      ## calculate ll null
      pl = plen[pki]
      ll_null[pki] <-  .log_loss_frag_ME(
        est_inser = lamb_all[pki],
        peak_length = pl,
        capturing_rates = capturing_rates,
        min_frag_length = min_frag_length,
        max_frag_length = max_frag_length,
        obs_pic_vec = pic_mat[pki,]
      )
      ll_full[pki] <- .log_loss_frag_ME(
        est_inser = lamb_full_1[pki],
        peak_length = pl,
        capturing_rates = capturing_rates_1,
        min_frag_length = min_frag_length,
        max_frag_length = max_frag_length,
        obs_pic_vec = pic_mat_1[pki,]
      ) +
        .log_loss_frag_ME(
          est_inser = lamb_full_2[pki],
          peak_length = pl,
          capturing_rates = capturing_rates_2,
          min_frag_length = min_frag_length,
          max_frag_length = max_frag_length,
          obs_pic_vec = pic_mat_2[pki,]
        )
    }
    p_val <- pchisq(2*( ll_full - ll_null), df = 1,lower.tail = F)

  }

  return(p_val)
}
