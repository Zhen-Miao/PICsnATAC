#' Get theoretical snATAC-seq distribution under condition 1
#'
#' @param insertion_rate The insertion rate (per 1000 base pairs)
#' @param peak_length The width of peak
#'
#' @return A vector of length 6, representing the probability of observing
#'  0 to >=5 counts
#'
#' @keywords internal
#' @noRd
.get_theoretical_c1 <- function(insertion_rate,
                                peak_length = 500) {
  lambda_ <- insertion_rate * peak_length / 1000
  p_W_m <- vector(mode = "numeric", length = 6)
  p_W_m[1] <- exp(-0.5 * lambda_) * (2 - exp(-0.5 * lambda_)) ## prob of 0
  p_W_m[2] <- exp(-0.5 * lambda_) * (lambda_ - 2 + 2 * exp(-0.5 * lambda_))
  p_W_m[3] <- exp(-1 * lambda_) *
    (exp(0.5 * lambda_) * (lambda_^2 - 4 * lambda_ + 8) - 8) / 4
  p_W_m[4] <- exp(-1 * lambda_) *
    (exp(0.5 * lambda_) *
       (lambda_^3 - 6 * lambda_^2 + 24 * lambda_ - 48) + 48) / 24
  p_W_m[5] <- exp(-1 * lambda_) *
    (exp(0.5 * lambda_) * (lambda_^4 - 8 * lambda_^3 + 48 *
      lambda_^2 - 192 * lambda_ + 384) + -384) / 192
  p_W_m[6] <- 1 - sum(p_W_m) ## prob of observing 5 or higher

  ## bi-allelic scenario
  p_W_m_conv <- stats::convolve(p_W_m, rev(p_W_m), type = "open")
  p_W_m_conv[6] <- sum(p_W_m_conv[6:11])

  ## return a vector of length 6, representing the probability of
  ## observing 0 to >=5 counts
  return(p_W_m_conv[1:6])
}


#' Construct an insertion-rate by mean count matrix
#'
#' @param insertion_rates A vector of insertion rates (per 1000 base pairs)
#' @param peak_lengths A vector containing width of peaks
#'
#' @return An insertion-rate by mean count matrix
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
      theo_mean_c1[iii, jjj] <- sum(p_W_s_theo * 0:5)
    }
  }
  rownames(theo_mean_c1) <- insertion_rates
  colnames(theo_mean_c1) <- peak_lengths

  return(theo_mean_c1)
}

#' Get theoretical distribution based on condition 1 and 2
#'
#' @param insertion_rate A single value of insertion rate
#' @param peak_length A single value of peak length, default 1000
#' @param pad_length Length of pad region, i.e., flanking regions of the peak
#'  where insertion can also happen
#' @param dirs Directions of insertion, for standard ATAC experiment (default),
#'  it is c(-1,1)
#' @param count_type output count type, default is 'PIC'
#' @param min_frag_length minimum fragment length, default 25
#' @param max_frag_length maximum fragment length, default 600
#' @param cap_insertion number of insertions capped at this value, default 20
#'
#' @return A vector of probability corresponding to different number of fragments
#' @noRd
#'
.get_theoretical_c12 <- function(
  insertion_rate,
  peak_length = 1000,
  pad_length = 0,
  dirs = c(-1, 1),
  count_type = "PIC",
  min_frag_length = 25,
  max_frag_length = 600,
  cap_insertion = 20 ## assume there will not be more than 20 insertions in
    ## one peak region
) {
  lambda_ <- insertion_rate * peak_length / 1000
  cap_fragment <- floor(cap_insertion / 2)
  p_W_m <- vector(length = cap_fragment)

  for (mi in seq_along(p_W_m)) {
    all_index_n <- (mi + 1):cap_insertion
    p_W_m_n <- vector(length = length(all_index_n))
    for (ni in seq_along(all_index_n)) {
      index_n <- all_index_n[ni]
      p_W_m_n_ni <- (choose(index_n - 1, mi) * 0.5^(index_n - 1) *
        lambda_^index_n * exp(-1 * lambda_)) /
        (factorial(index_n))
      if (length(p_W_m_n_ni) != 1L || !is.finite(p_W_m_n_ni)) {
        stop(
          "Could not compute a finite theoretical probability for insertion rate ",
          insertion_rate,
          call. = FALSE
        )
      }
      p_W_m_n[ni] <- p_W_m_n_ni
    }
    p_W_m[mi] <- sum(p_W_m_n)
  }


  p_W_s <- vector(length = cap_fragment)

  exp_capture <- (exp(
    -1 * min_frag_length * insertion_rate / 1000
  ) -
    exp(-1 * max_frag_length * insertion_rate / 1000)
  ) /
    (1 - exp(-1 * peak_length * insertion_rate / 1000))

  for (ki in seq_along(p_W_s)) {
    all_index_m <- ki:cap_fragment
    p_W_s_k <- vector(length = length(all_index_m))
    for (mii in seq_along(p_W_s_k)) {
      index_m <- all_index_m[mii]
      p_W_s_k[mii] <- choose(index_m, ki) * exp_capture^ki *
        (1 - exp_capture)^(index_m - ki) * p_W_m[index_m]
    }
    p_W_s[ki] <- sum(p_W_s_k)
  }

  ## bi-allelic scenario, we convolve it with itself
  p_W_s <- c(1 - sum(p_W_s), p_W_s)
  p_W_s_conv <- stats::convolve(p_W_s, rev(p_W_s), type = "open")
  return(p_W_s_conv[2:(cap_fragment + 1)])
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
#' @return An insertion rate by peak length matrix, where each
#' element of the matrix corresponds to the number of expected fragments
#' @noRd
.insertion_to_c12 <- function(
    insertion_rates = (1:2000) * 0.01,
    peak_lengths = 4:20 * 50,
    pad_length = 0,
    dirs = c(-1, 1),
    count_type = "fragment",
    min_frag_length = 25,
    max_frag_length = 600,
    cap_insertion = 20) {
  ## initialize the matrix
  theo_mean_c12 <- matrix(
    nrow = length(insertion_rates),
    ncol = length(peak_lengths)
  )

  for (iii in seq_along(insertion_rates)) {
    for (jjj in seq_along(peak_lengths)) {
      p_W_s_theo <- .get_theoretical_c12(
        insertion_rate = insertion_rates[iii],
        peak_length = peak_lengths[jjj],
        pad_length = pad_length,
        dirs = dirs,
        count_type = count_type,
        min_frag_length = min_frag_length,
        max_frag_length = max_frag_length,
        cap_insertion = cap_insertion
      )
      theo_mean_c12[iii, jjj] <- sum(p_W_s_theo * seq_along(p_W_s_theo))
    }
  }
  rownames(theo_mean_c12) <- insertion_rates
  colnames(theo_mean_c12) <- peak_lengths

  return(theo_mean_c12)
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
#'
#' @keywords internal
#' @noRd
.log_loss_frag_c1 <- function(est_inser,
                              peak_length,
                              cap_fragment = 5,
                              capturing_rates,
                              obs_pic_vec) {
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
        pt[as.character(ot), as.character(k)] <-
          choose(k, ot) * ct^ot * (1 - ct)^(k - ot) * p_W_s_theo[k + 1]
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


#' Validate a peak-by-cell count matrix
#'
#' @param pic_mat A base or sparse matrix.
#'
#' @return `pic_mat`, invisibly.
#' @noRd
.validate_pic_matrix <- function(pic_mat) {
  if (!is.matrix(pic_mat) && !methods::is(pic_mat, "Matrix")) {
    stop("pic_mat must be a matrix or Matrix object", call. = FALSE)
  }
  if (nrow(pic_mat) == 0L || ncol(pic_mat) == 0L) {
    stop("pic_mat must have at least one peak and one cell", call. = FALSE)
  }
  if (anyNA(pic_mat) || any(!is.finite(pic_mat)) || any(pic_mat < 0)) {
    stop("pic_mat must contain finite, non-negative counts", call. = FALSE)
  }
  invisible(pic_mat)
}


#' Discretize capture rates onto the likelihood model grid
#'
#' @param capturing_rates A numeric vector of probabilities.
#' @param n_cells Expected vector length.
#'
#' @return A numeric vector on the grid used by `.log_loss_frag_c1()`.
#' @noRd
.discretize_capturing_rates <- function(capturing_rates, n_cells) {
  if (!is.numeric(capturing_rates) || length(capturing_rates) != n_cells ||
      anyNA(capturing_rates) || any(!is.finite(capturing_rates)) ||
      any(capturing_rates < 0 | capturing_rates > 1)) {
    stop(
      "capturing_rates must contain one finite probability in [0, 1] per cell",
      call. = FALSE
    )
  }

  capturing_rates <- ceiling(capturing_rates * 10) / 10
  capturing_rates[capturing_rates < 0.2] <- 0.2
  capturing_rates[capturing_rates > 0.9] <- 0.99
  capturing_rates
}


#' Validate the requested parallelism
#'
#' @param n_cores Number of worker cores.
#'
#' @return An integer scalar.
#' @noRd
.validate_n_cores <- function(n_cores) {
  if (length(n_cores) != 1L || !is.numeric(n_cores) || !is.finite(n_cores) ||
      n_cores < 1L || n_cores != as.integer(n_cores)) {
    stop("n_cores must be a positive integer", call. = FALSE)
  }
  as.integer(n_cores)
}


#' Apply a function across peaks with a cross-platform fallback
#'
#' @param x Values to iterate over.
#' @param fun Function to apply.
#' @param n_cores Number of requested cores.
#'
#' @return A list.
#' @noRd
.parallel_peak_lapply <- function(x, fun, n_cores) {
  n_cores <- .validate_n_cores(n_cores)
  if (n_cores == 1L) {
    return(lapply(x, fun))
  }
  if (.Platform$OS.type == "windows") {
    warning(
      "Fork-based parallelism is unavailable on Windows; using one core",
      call. = FALSE
    )
    return(lapply(x, fun))
  }
  parallel::mclapply(x, fun, mc.cores = n_cores)
}


#' Parse peak widths from matrix row names
#'
#' @param peak_names Names in `seqname:start-end` form.
#'
#' @return A numeric vector of peak widths.
#' @noRd
.peak_widths_from_names <- function(peak_names) {
  if (is.null(peak_names) || anyNA(peak_names) || any(!nzchar(peak_names))) {
    stop(
      "pic_mat must have non-missing seqname:start-end row names when plen is NULL",
      call. = FALSE
    )
  }

  parsed <- utils::strcapture(
    "^.+:([0-9]+)-([0-9]+)$",
    peak_names,
    proto = list(start = numeric(), end = numeric())
  )
  if (anyNA(parsed$start) || anyNA(parsed$end) ||
      any(parsed$end <= parsed$start)) {
    stop(
      "pic_mat row names must use seqname:start-end with end greater than start",
      call. = FALSE
    )
  }
  parsed$end - parsed$start
}


#' Apply the DAR count-quality policy
#'
#' @param pic_mat A base or sparse peak-by-cell count matrix.
#' @param count_cap Largest count represented explicitly by the model.
#' @param artifact_threshold Counts strictly above this heuristic threshold are
#'  treated as unreliable and replaced by zero. Use `Inf` to disable this rule.
#'
#' @return A matrix of the same broad representation (base or sparse), with
#'  artifacts replaced by zero and remaining high counts capped.
#' @noRd
.prepare_pic_counts <- function(pic_mat,
                                count_cap = 5L,
                                artifact_threshold = 10L) {
  .validate_pic_matrix(pic_mat)
  if (length(count_cap) != 1L || !is.numeric(count_cap) ||
      !is.finite(count_cap) || count_cap < 1L ||
      count_cap != as.integer(count_cap)) {
    stop("count_cap must be a positive integer", call. = FALSE)
  }
  if (length(artifact_threshold) != 1L ||
      !is.numeric(artifact_threshold) || is.na(artifact_threshold) ||
      artifact_threshold < count_cap) {
    stop(
      "artifact_threshold must be at least count_cap (or Inf to disable it)",
      call. = FALSE
    )
  }

  is_sparse <- methods::is(pic_mat, "Matrix")
  if (is_sparse) {
    pic_mat <- methods::as(Matrix::Matrix(pic_mat, sparse = TRUE), "dMatrix")
    counts <- pic_mat@x
  } else {
    counts <- as.numeric(pic_mat)
  }
  if (any(counts != floor(counts))) {
    stop("pic_mat counts must be whole numbers", call. = FALSE)
  }

  is_artifact <- counts > artifact_threshold
  counts[is_artifact] <- 0
  counts[!is_artifact & counts > count_cap] <- count_cap

  if (is_sparse) {
    pic_mat@x <- counts
    return(Matrix::drop0(pic_mat))
  }

  matrix(
    counts,
    nrow = nrow(pic_mat),
    ncol = ncol(pic_mat),
    dimnames = dimnames(pic_mat)
  )
}


#' Calculate maximized log-likelihoods from observed PIC counts
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param plen A vector of peak widths
#' @param n_cores A positive integer specifying the number of cores. On Windows,
#'  values greater than one fall back to serial evaluation.
#'
#' @return A numeric vector containing the maximized log-likelihood for each
#'  peak over insertion rates from 0.01 to 20.
#' @export
#'
obs_to_insertion_MLE_obj <- function(pic_mat,
                                     capturing_rates,
                                     plen,
                                     n_cores = 1L) {
  .validate_pic_matrix(pic_mat)
  n_cores <- .validate_n_cores(n_cores)
  if (!is.numeric(plen) || length(plen) != nrow(pic_mat) ||
      anyNA(plen) || any(!is.finite(plen)) || any(plen <= 0)) {
    stop("plen must contain one positive, finite width per peak", call. = FALSE)
  }
  capturing_rates <- .discretize_capturing_rates(
    capturing_rates,
    ncol(pic_mat)
  )
  observed_counts <- if (methods::is(pic_mat, "Matrix")) {
    numeric_pic_mat <- methods::as(
      Matrix::Matrix(pic_mat, sparse = TRUE),
      "dMatrix"
    )
    numeric_pic_mat@x
  } else {
    pic_mat
  }
  if (any(observed_counts != floor(observed_counts)) ||
      any(observed_counts > 5)) {
    stop("pic_mat must contain whole-number counts between 0 and 5",
      call. = FALSE
    )
  }

  n_para <- length(plen)

  ## iterations
  optim_results <- .parallel_peak_lapply(seq_len(n_para), function(pp) {
    stats::optimize(
      f = .log_loss_frag_c1,
      interval = c(0.01, 20),
      peak_length = plen[pp],
      capturing_rates = capturing_rates,
      obs_pic_vec = as.numeric(pic_mat[pp, , drop = TRUE]),
      maximum = TRUE
    )$objective
  }, n_cores = n_cores)

  unlist(optim_results, use.names = FALSE)
}


#' Moment estimator for insertion rates from observed values
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param cell_type_labels A vector of cell type labels for each cell
#' @param cap_insertion The maximum number of insertions in a peak region.
#' @param insertion_rates The range of insertion rates (per 1000 bp) to be
#'  considered.
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
    insertion_rates = (1:2000) * 0.01,
    peak_lengths = 4:20 * 50) {
  .validate_pic_matrix(pic_mat)
  if (length(capturing_rates) != ncol(pic_mat) ||
      !is.numeric(capturing_rates) || anyNA(capturing_rates) ||
      any(!is.finite(capturing_rates)) ||
      any(capturing_rates < 0 | capturing_rates > 1)) {
    stop(
      "capturing_rates must contain one finite probability in [0, 1] per cell",
      call. = FALSE
    )
  }
  if (length(cell_type_labels) != ncol(pic_mat) || anyNA(cell_type_labels)) {
    stop("cell_type_labels must contain one non-missing label per cell",
      call. = FALSE
    )
  }
  cell_type_labels <- as.character(cell_type_labels)
  if (!is.numeric(insertion_rates) || length(insertion_rates) == 0L ||
      anyNA(insertion_rates) || any(!is.finite(insertion_rates)) ||
      any(insertion_rates <= 0)) {
    stop("insertion_rates must contain positive, finite values", call. = FALSE)
  }
  if (!is.numeric(peak_lengths) || length(peak_lengths) < 2L ||
      anyNA(peak_lengths) || any(!is.finite(peak_lengths)) ||
      any(diff(peak_lengths) <= 0) ||
      !isTRUE(all.equal(diff(peak_lengths), rep(diff(peak_lengths)[1L],
        length(peak_lengths) - 1L)))) {
    stop("peak_lengths must be an equally spaced, increasing numeric vector",
      call. = FALSE
    )
  }
  if (length(cap_insertion) != 1L || !is.numeric(cap_insertion) ||
      !is.finite(cap_insertion) || cap_insertion < 2L ||
      cap_insertion != as.integer(cap_insertion)) {
    stop("cap_insertion must be an integer of at least 2", call. = FALSE)
  }
  if (length(min_frag_length) != 1L || length(max_frag_length) != 1L ||
      !is.finite(min_frag_length) || !is.finite(max_frag_length) ||
      min_frag_length < 0 || max_frag_length <= min_frag_length) {
    stop(
      "min_frag_length and max_frag_length must define a positive interval",
      call. = FALSE
    )
  }

  ## get the matrix
  c12_mat <- .insertion_to_c12(
    insertion_rates = insertion_rates,
    peak_lengths = peak_lengths,
    min_frag_length = min_frag_length,
    max_frag_length = max_frag_length,
    cap_insertion = cap_insertion
  )

  ## Group peaks by length.
  plen <- .peak_widths_from_names(rownames(pic_mat))

  ## assign group
  plen_gap <- peak_lengths[2] - peak_lengths[1]
  plen_group <- ceiling((plen - peak_lengths[1]) / plen_gap) + 1
  plen_group[plen_group < 1L] <- 1L
  plen_group[plen_group > length(peak_lengths)] <- length(peak_lengths)

  ## to prevent the low acc cells, adjust capturing rate
  capturing_rates[capturing_rates < 0.2] <- 0.2

  ## do this for each cell type
  cell_types <- unique(cell_type_labels)
  plen_group <- as.character(plen_group)
  pgroups <- unique(plen_group)

  est_inser_all_matrix <- matrix(
    NA_real_,
    nrow = nrow(pic_mat),
    ncol = length(cell_types),
    dimnames = list(rownames(pic_mat), cell_types)
  )

  for (ct_index in seq_along(cell_types)) {
    ct <- cell_types[ct_index]
    sel_cell <- cell_type_labels == ct
    pic_sub <- pic_mat[, sel_cell, drop = FALSE]
    w_bar <- Matrix::rowSums(
      pic_sub %*% Matrix::Diagonal(x = 1 / capturing_rates[sel_cell])
    )
    w_bar <- w_bar / sum(sel_cell)

    for (pg in pgroups) {
      c12_sub <- c12_mat[, as.numeric(pg)]
      peak_selection <- plen_group == pg
      w_bar_sub <- w_bar[peak_selection]
      lkup <- outer(w_bar_sub, c12_sub, "-")
      lkup <- abs(lkup)

      # for each row, look up the smallest value
      est_inser_ind <- apply(lkup, 1, which.min)
      est_inser_all_matrix[peak_selection, ct_index] <-
        insertion_rates[est_inser_ind]
    }
  }

  est_inser_all_matrix
}


#' Log loss from fragment count
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
    obs_pic_vec) {
  n_cells <- length(obs_pic_vec)

  p_W_s_theo <- .get_theoretical_c12(
    insertion_rate = est_inser,
    peak_length = peak_length,
    min_frag_length = min_frag_length,
    max_frag_length = max_frag_length
  )

  ## include p(W_s = 0)
  p_W_s_0 <- 1 - sum(p_W_s_theo)
  p_W_s_theo <- c(p_W_s_0, p_W_s_theo)
  lg_p_W_o_t <- vector(mode = "numeric", length = n_cells)
  for (cell_i in seq_len(n_cells)) {
    ot <- obs_pic_vec[cell_i]
    ct <- capturing_rates[cell_i]
    pt <- vector(mode = "numeric", length = cap_fragment - ot + 1)
    names(pt) <- ot:cap_fragment
    for (k in ot:cap_fragment) {
      pt[as.character(k)] <- choose(k, ot) *
        ct^ot * (1 - ct)^(k - ot) * p_W_s_theo[k + 1]
    }
    lg_p_W_o_t[cell_i] <- log(sum(pt))
  }
  ## remove the effect of very large value
  lg_p_W_o_t[lg_p_W_o_t < -15] <- -15
  ll <- sum(lg_p_W_o_t)
  return(ll)
}



#' Calculate insertion rates from observed PIC counts
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param plen A vector of peak widths
#' @param n_cores A positive integer specifying the number of cores. On Windows,
#'  values greater than one fall back to serial evaluation.
#'
#' @return The optimized insertion rate (per 1000 base pairs) over
#'  insertion rates from 0.01 to 20
#' @export
#'
obs_to_insertion_MLE_lam <- function(pic_mat,
                                     capturing_rates,
                                     plen,
                                     n_cores = 1L) {
  .validate_pic_matrix(pic_mat)
  n_cores <- .validate_n_cores(n_cores)
  if (!is.numeric(plen) || length(plen) != nrow(pic_mat) ||
      anyNA(plen) || any(!is.finite(plen)) || any(plen <= 0)) {
    stop("plen must contain one positive, finite width per peak", call. = FALSE)
  }
  capturing_rates <- .discretize_capturing_rates(
    capturing_rates,
    ncol(pic_mat)
  )
  observed_counts <- if (methods::is(pic_mat, "Matrix")) {
    numeric_pic_mat <- methods::as(
      Matrix::Matrix(pic_mat, sparse = TRUE),
      "dMatrix"
    )
    numeric_pic_mat@x
  } else {
    pic_mat
  }
  if (any(observed_counts != floor(observed_counts)) ||
      any(observed_counts > 5)) {
    stop("pic_mat must contain whole-number counts between 0 and 5",
      call. = FALSE
    )
  }

  n_para <- length(plen)

  ## iterations
  optim_results <- .parallel_peak_lapply(seq_len(n_para), function(pp) {
    stats::optimize(
      f = .log_loss_frag_c1,
      interval = c(0.01, 20),
      peak_length = plen[pp],
      capturing_rates = capturing_rates,
      obs_pic_vec = as.numeric(pic_mat[pp, , drop = TRUE]),
      maximum = TRUE
    )$maximum
  }, n_cores = n_cores)
  unlist(optim_results, use.names = FALSE)
}

#' Compute p value for DAR test between two cell types
#'
#' @param pic_mat The observed peak by cell PIC count matrix
#' @param capturing_rates A vector of estimated capturing rates for each cell
#' @param cell_type_labels A vector specifying cell type labels
#' @param estimation_approach The approach for parameter estimation,
#'   either 'MLE'
#'   for condition 1 or 'ME' for condition 1+2.
#'   The 'MLE' approach is more accurate
#'   and usually it has a higher power, but it ignores the size filtering step
#'   in snATAC-seq data generation. Default is 'MLE'
#' @param n_cores A positive integer specifying the number of cores. On Windows,
#'  values greater than one fall back to serial evaluation. Default = 1.
#' @param plen A vector of peak length
#' @param min_frag_length The value for the s1 hyperparameter in the ssPoisson
#'   distribution, this stands for the minimum fragment length requirement such
#'   that the fragment can be amplifiable and mappable to genome. Default = 25
#' @param max_frag_length The value for the s2 hyperparameter in the ssPoisson
#'   distribution, this stands for the max fragment length requirement such
#'   that the fragment can be amplifiable. Default = 600
#' @param artifact_threshold A heuristic upper threshold for unreliable counts.
#'   Counts strictly greater than this value may reflect mapping errors or
#'   unusual fragment structures and are replaced by zero before the remaining
#'   counts are capped at 5 for the likelihood model. The default is 10; use
#'   `Inf` to disable artifact replacement. Because this cutoff is assay- and
#'   pipeline-dependent, sensitivity analyses with alternative values are
#'   recommended.
#'
#' @return A numeric vector containing one likelihood-ratio-test p-value per
#'   peak, named from `rownames(pic_mat)` when available.
#' @export
#'
DAR_by_LRT <- function(pic_mat,
                       capturing_rates,
                       cell_type_labels,
                       n_cores = 1L,
                       plen = NULL,
                       min_frag_length = 25,
                       max_frag_length = 600,
                       estimation_approach = "MLE",
                       artifact_threshold = 10L) {
  .validate_pic_matrix(pic_mat)
  n_cores <- .validate_n_cores(n_cores)
  if (!is.character(estimation_approach) ||
      length(estimation_approach) != 1L ||
      !(estimation_approach %in% c("MLE", "ME"))) {
    stop("estimation_approach must be either 'MLE' or 'ME'", call. = FALSE)
  }
  if (length(cell_type_labels) != ncol(pic_mat) || anyNA(cell_type_labels)) {
    stop("cell_type_labels must contain one non-missing label per cell",
      call. = FALSE
    )
  }
  cell_type_labels <- as.character(cell_type_labels)
  ct_uniq <- unique(cell_type_labels)
  if (length(ct_uniq) != 2L) {
    stop("DAR_by_LRT requires exactly two cell types", call. = FALSE)
  }
  capturing_rates <- .discretize_capturing_rates(
    capturing_rates,
    ncol(pic_mat)
  )

  if (as.double(nrow(pic_mat)) * as.double(ncol(pic_mat)) >= (2^31 - 1)) {
    stop(
      "pic_mat is too large for a single test; split it by peaks and combine the results",
      call. = FALSE
    )
  }

  ## save some values
  n_pks <- nrow(pic_mat)

  ## pk length
  if (is.null(plen)) {
    plen <- .peak_widths_from_names(rownames(pic_mat))
    ## we only need to group peak length by n*100 bp
    plen <- ceiling(plen / 100) * 100
  } else if (!is.numeric(plen) || length(plen) != n_pks ||
      anyNA(plen) || any(!is.finite(plen)) || any(plen <= 0)) {
    stop("plen must contain one positive, finite width per peak", call. = FALSE)
  }

  ## Counts above the model support are capped. Counts strictly above the
  ## configurable quality threshold are first treated as unreliable artifacts.
  pic_mat <- .prepare_pic_counts(
    pic_mat,
    count_cap = 5L,
    artifact_threshold = artifact_threshold
  )

  ## split the matrix by cell types
  pic_mat_1 <- pic_mat[, cell_type_labels == ct_uniq[1], drop = FALSE]
  pic_mat_2 <- pic_mat[, cell_type_labels == ct_uniq[2], drop = FALSE]

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
      capturing_rates = capturing_rates_1,
      plen = plen,
      n_cores = n_cores
    )

    ll_full_2_mle <- obs_to_insertion_MLE_obj(
      pic_mat = pic_mat_2,
      capturing_rates = capturing_rates_2,
      plen = plen,
      n_cores = n_cores
    )

    ## p value is obtained by chi-squared statistics
    p_val <- stats::pchisq(2 * (ll_full_1_mle + ll_full_2_mle - ll_all_mle),
      df = 1, lower.tail = FALSE
    )
  } else if (estimation_approach == "ME") {
    ## calculate estimated lambda for null hypothesis
    lamb_all <- obs_to_insertion_ME(
      pic_mat = pic_mat,
      capturing_rates = capturing_rates,
      min_frag_length = min_frag_length,
      max_frag_length = max_frag_length,
      cell_type_labels = rep("A", length = length(capturing_rates))
    )
    lamb_all <- lamb_all[, 1]

    if (anyNA(lamb_all)) {
      stop("Moment estimation produced missing values under the null model",
        call. = FALSE
      )
    }

    ## calculate estimated lambda for alternative hypothesis
    lamb_full <- obs_to_insertion_ME(
      pic_mat = pic_mat,
      capturing_rates = capturing_rates,
      min_frag_length = min_frag_length,
      max_frag_length = max_frag_length,
      cell_type_labels = cell_type_labels
    )
    lamb_full_1 <- lamb_full[, 1]
    lamb_full_2 <- lamb_full[, 2]


    if (anyNA(lamb_full_1) || anyNA(lamb_full_2)) {
      stop("Moment estimation produced missing values under the full model",
        call. = FALSE
      )
    }

    ## test for each peak
    ll_null <- vector(mode = "numeric", length = n_pks)
    p_val <- ll_full <- ll_null


    for (pki in seq_len(n_pks)) {
      ## calculate ll null
      pl <- plen[pki]
      ll_null[pki] <- .log_loss_frag_ME(
        est_inser = lamb_all[pki],
        peak_length = pl,
        capturing_rates = capturing_rates,
        min_frag_length = min_frag_length,
        max_frag_length = max_frag_length,
        obs_pic_vec = as.numeric(pic_mat[pki, , drop = TRUE])
      )
      ll_full[pki] <- .log_loss_frag_ME(
        est_inser = lamb_full_1[pki],
        peak_length = pl,
        capturing_rates = capturing_rates_1,
        min_frag_length = min_frag_length,
        max_frag_length = max_frag_length,
        obs_pic_vec = as.numeric(pic_mat_1[pki, , drop = TRUE])
      ) +
        .log_loss_frag_ME(
          est_inser = lamb_full_2[pki],
          peak_length = pl,
          capturing_rates = capturing_rates_2,
          min_frag_length = min_frag_length,
          max_frag_length = max_frag_length,
          obs_pic_vec = as.numeric(pic_mat_2[pki, , drop = TRUE])
        )
    }
    p_val <- stats::pchisq(
      2 * (ll_full - ll_null),
      df = 1,
      lower.tail = FALSE
    )
  }

  names(p_val) <- rownames(pic_mat)
  p_val
}
