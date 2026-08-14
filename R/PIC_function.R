## PIC function


#' Function for loading fragment files and filter by cell barcodes
#'
#' @param fragment_tsv_gz_file_location The 10X Cell Ranger output
#'  fragment.tsv.gz file location. This can usually be found at the /out
#'  directory from Cell Ranger output
#' @param cells Cell barcode labels as a character vector.
#' @param verbose Whether to output progress message. Default TRUE
#'
#' @return A data frame containing fragments filtered by cell barcode.
#' @export
load_fragments <- function(
    fragment_tsv_gz_file_location, cells, verbose = TRUE) {
  if (!is.character(fragment_tsv_gz_file_location) ||
      length(fragment_tsv_gz_file_location) != 1L ||
      is.na(fragment_tsv_gz_file_location) ||
      !file.exists(fragment_tsv_gz_file_location)) {
    stop("fragment_tsv_gz_file_location must identify an existing file",
      call. = FALSE
    )
  }
  if (!is.character(cells) || length(cells) == 0L || anyNA(cells) ||
      any(!nzchar(cells))) {
    stop("cells must be a non-empty character vector without missing values",
      call. = FALSE
    )
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }

  f1 <- data.table::fread(fragment_tsv_gz_file_location,
    header = FALSE,
    select = 1:4
  )
  ## data.table show inconsistent performance
  # setnames(f1, c('seqname','start','end','cell_barcode'))
  # f1 <- f1[f1$cell_barcode %in% cells]

  ## convert to data.frame format
  f1 <- as.data.frame(f1)
  colnames(f1) <- c("seqname", "start", "end", "cell_barcode")

  cells_retain <- f1$cell_barcode %in% cells
  n_fragments_retained <- sum(cells_retain)

  ## report the proportion of reads in cell barcodes
  prop_bcreads <- n_fragments_retained / nrow(f1)
  if (verbose) {
    cat(sprintf("proportion of reads in cell barcodes is %.2f\n", prop_bcreads))
  }

  ## error when no cells found in the fragment file
  if (n_fragments_retained < 1L) {
    stop("Cell barcodes not found in fragment file; please check the input",
      call. = FALSE
    )
  } else if (n_fragments_retained < 10L) {
    warning(
      "Fewer than 10 fragments matched the requested cell barcodes; check the input",
      call. = FALSE
    )
  }

  f1 <- f1[cells_retain, ]
  f1
}


#' Count number of paired insertions in each peak
#'
#' @param peak_sets GRanges object of peak sets
#' @param filtered_fragments filtered fragment also as a GRanges object
#' @param n_features Number of features (peaks)
#' @param extend_size How long should we extend the exact insertion
#' site as accessible window
#'
#' @return A sparse vector of PIC for each peak
#' @export
count_peaks <- function(peak_sets, filtered_fragments,
                        extend_size, n_features = length(peak_sets)) {
  if (length(n_features) != 1L || is.na(n_features) ||
      n_features != length(peak_sets)) {
    stop("n_features must equal length(peak_sets)", call. = FALSE)
  }

  ## get start and end position
  f1s <- GenomicRanges::resize(filtered_fragments, width = 1, fix = "start")
  f1e <- GenomicRanges::resize(filtered_fragments, width = 1, fix = "end")

  ## extend the exact insertion site
  f1s <- GenomicRanges::resize(f1s, width = extend_size, fix = "center")
  f1e <- GenomicRanges::resize(f1e, width = extend_size, fix = "center")

  overlaped_s <- GenomicRanges::findOverlaps(f1s, peak_sets, select = "first")
  overlaped_e <- GenomicRanges::findOverlaps(f1e, peak_sets, select = "first")
  overlaped_e_nodc <- overlaped_e

  overlaped_e_nodc[overlaped_s == overlaped_e_nodc] <- NA
  ol <- c(overlaped_s, overlaped_e_nodc)

  counts <- as.data.frame(table(ol), stringsAsFactors = FALSE)
  counts$ol <- as.integer(counts$ol)

  out_vec <- Matrix::sparseVector(
    x = counts$Freq,
    i = counts$ol,
    length = n_features
  )

  out_vec
}


#' Generate a sparse matrix from a list of sparse vectors
#'
#' @param chunk A chunk of lists with each element being a sparseVector
#'
#' @return A sparseMatrix that is a column bind of all sparseVectors within
#'  the chunk
#' @noRd
make_s_mat_from_s_vec <- function(chunk, n_features) {
  indices <- lapply(chunk, function(y) y@i)
  values <- lapply(chunk, function(y) y@x)

  i <- unlist(indices, use.names = FALSE)
  j <- unlist(lapply(
    seq_along(chunk),
    function(k) rep(k, length(indices[[k]]))
  ), use.names = FALSE)
  x <- unlist(values, use.names = FALSE)

  s_mat <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(n_features, length(chunk))
  )

  ## add cell names to the matrix
  colnames(s_mat) <- names(chunk)
  s_mat
}


#' From a (potentially large) list of sparseVectors into a full sparseMatrix
#'
#' @param list_s_vetors A chunk of lists with each element being a sparseVector
#' @param n_features Number of features (peaks)
#'
#'
#' @return A sparseMatrix that is a column bind of all sparseVectors
#' @noRd
list_to_sparseMatrix <- function(list_s_vetors, n_features) {
  n_cells <- as.numeric(length(list_s_vetors))
  if (n_cells < 1L) {
    stop("list_s_vetors must contain at least one cell", call. = FALSE)
  }
  if (length(n_features) != 1L || !is.finite(n_features) || n_features < 0) {
    stop("n_features must be a non-negative scalar", call. = FALSE)
  }

  max_cells_per_chunk <- if (n_features == 0L) {
    n_cells
  } else {
    max(1, floor((2^31 - 1) / n_features))
  }

  if (n_cells > max_cells_per_chunk) {
    chunks <- split(
      list_s_vetors,
      ceiling(seq_along(list_s_vetors) / max_cells_per_chunk)
    )
    sparse_matrices <- lapply(
      chunks,
      make_s_mat_from_s_vec,
      n_features = n_features
    )
    sparse_matrices <- do.call(cbind, sparse_matrices)
  } else {
    sparse_matrices <- make_s_mat_from_s_vec(list_s_vetors, n_features)
  }
  sparse_matrices
}


#' convert peak_set into GRanges object if input is a data.frame
#'
#' @param peak_sets A data.frame object of peaks that we want to use as features
#'  The first column should be seqname, (e.g., chr1); the second column should
#'  be the start site, and the third column should be the end site. This should
#'  be after 5 bp and 4 bp correction of Tn5 insertion location.
#'
#' @return A GRanges object that contain the same information as peak_sets
#' @export
data_frame_to_GRanges <- function(peak_sets) {
  if (!is.data.frame(peak_sets) && !is.matrix(peak_sets)) {
    stop("peak_sets must be a data.frame or matrix", call. = FALSE)
  }
  if (ncol(peak_sets) < 3L) {
    stop("peak_sets must contain at least three columns", call. = FALSE)
  }

  peak_sets <- as.data.frame(peak_sets)
  ## if colnames not specified, we specify by order
  if (is.null(colnames(peak_sets)) ||
    !(all(c("seqname", "start", "end") %in% colnames(peak_sets)))) {
    colnames(peak_sets)[1:3] <- c("seqname", "start", "end")
  }

  ## convert into GRanges
  peak_sets <- tryCatch(
    GenomicRanges::makeGRangesFromDataFrame(peak_sets),
    error = function(e) {
      stop(
        "Could not convert peak_sets to a GRanges object: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  peak_sets
}


#' Format peak coordinates for matrix row names
#'
#' @param peak_sets A `GRanges` object.
#'
#' @return A character vector in `seqname:start-end` form.
#' @noRd
.format_peak_names <- function(peak_sets) {
  paste0(
    as.character(GenomeInfoDb::seqnames(peak_sets)),
    ":",
    BiocGenerics::start(peak_sets),
    "-",
    BiocGenerics::end(peak_sets)
  )
}



#' Count snATAC-seq data matrix with Paired Insertion Counting (PIC)
#'
#' @param cells Cell barcode labels as a character vector.
#' @param fragment_tsv_gz_file_location The 10X Cell Ranger output
#'  fragment.tsv.gz file location. This can usually be found at the /out
#'  directory from Cell Ranger output
#' @param peak_sets The set of peaks as a GenomicRanges object. This will be
#'  the features for the data matrix. Alternatively, this can be a data.frame
#'  and the function will convert it into a GenomicRanges object
#' @param deduplicate Whether to include deduplicate step where within
#'  the same cell,
#'  fragments with identical start and end location will be deduplicated.
#'  This is usually unnecessary for Cell Ranger ATAC output, since
#'  Cell Ranger ATAC has already deduplicated the fragments.
#'  But for dsc-ATAC-seq data, this step will
#'  be helpful and recommended.
#' @param load_full Whether to load the whole fragment.tsv.gz file into memory.
#'  If set to `FALSE`, the function loads it by chromosome to save RAM. This
#'  mode requires a block-gzipped file and its Tabix index (`.tbi`).
#' @param extend_size How long should we extend the exact insertion site as
#'  an accessible window, in base pairs.
#' @param verbose Whether to output progress information including the progress
#'  bar
#'
#' @return A sparse peak-by-cell PIC count matrix.
#' @export
#'
PIC_counting <- function(cells,
                         fragment_tsv_gz_file_location,
                         peak_sets,
                         deduplicate = FALSE,
                         load_full = TRUE,
                         extend_size = 5L,
                         verbose = TRUE) {
  ## check input
  if (!is.character(cells) || length(cells) == 0L || anyNA(cells) ||
      any(!nzchar(cells))) {
    stop("cells must be a non-empty character vector without missing values",
      call. = FALSE
    )
  }
  if (anyDuplicated(cells)) {
    stop("cells must not contain duplicate barcodes", call. = FALSE)
  }
  if (length(extend_size) != 1L || !is.finite(extend_size) ||
      extend_size < 0 || extend_size != as.integer(extend_size)) {
    stop("extend_size must be a non-negative integer", call. = FALSE)
  }
  if (!is.logical(deduplicate) || length(deduplicate) != 1L ||
      is.na(deduplicate)) {
    stop("deduplicate must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(load_full) || length(load_full) != 1L || is.na(load_full)) {
    stop("load_full must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.character(fragment_tsv_gz_file_location) ||
      length(fragment_tsv_gz_file_location) != 1L ||
      is.na(fragment_tsv_gz_file_location) ||
      !file.exists(fragment_tsv_gz_file_location)) {
    stop("fragment_tsv_gz_file_location must identify an existing file",
      call. = FALSE
    )
  }

  ## we accept peak_sets to be a GRanges or we convert it into one
  if (!methods::is(peak_sets, "GRanges")) {
    peak_sets <- data_frame_to_GRanges(peak_sets)
  }

  n_cells <- length(cells)
  n_features <- length(peak_sets)

  if (n_features == 0L) {
    stop("peak_sets must contain at least one peak", call. = FALSE)
  }

  ## if load full files
  if (load_full) {
    ## create output object
    out_summ <- rep(list(), length = n_cells)
    names(out_summ) <- cells

    ## load fragment files
    f1 <- load_fragments(
      fragment_tsv_gz_file_location = fragment_tsv_gz_file_location,
      cells = cells,
      verbose = verbose
    )

    ## require the end to be larger than start -- this is useful for
    ## s3-ATAC-seq data if end smaller than start
    if (sum(f1$"start" - 1 >= f1$"end") >= 1) {
      f1_sub1 <- f1[f1$"start" - 1 < f1$"end", ]
      f1_sub2 <- f1[f1$"start" - 1 >= f1$"end", ]
      f1_sub2s <- f1_sub2$"start"
      f1_sub2$"start" <- f1_sub2$"end"
      f1_sub2$"end" <- f1_sub2s
      f1 <- rbind(f1_sub1, f1_sub2)
      rm(f1_sub2)
      rm(f1_sub2s)
      rm(f1_sub1)
    }


    ## generate GenomicRanges object
    f1 <- GenomicRanges::makeGRangesFromDataFrame(f1,
      keep.extra.columns = TRUE
    )

    ## pre-sort fragments
    f1_s <- IRanges::subsetByOverlaps(f1,
      ranges = peak_sets,
      maxgap = ceiling(extend_size / 2)
    )
    rm(f1)

    n_subset <- ceiling(n_cells / 500)
    f1k <- rep(list(), length = n_subset)
    for (i in seq_len(n_subset)) {
      s <- (i - 1) * 500 + 1
      e <- min(i * 500, n_cells)

      ## deduplicate f1_s
      if (deduplicate) {
        f1k[[i]] <- unique(f1_s[f1_s$cell_barcode %in% cells[s:e], ])
      } else {
        f1k[[i]] <- f1_s[f1_s$cell_barcode %in% cells[s:e], ]
      }
    }
    rm(f1_s)
    gc()

    ## progress bar
    if (verbose) {
      pb <- progress::progress_bar$new(
        total = n_cells,
        format = "[:bar] :percent finished, elapsed: :elapsed",
        clear = FALSE,
        width = 60
      )
      cat("Computing peak vector for each cell.\n")
    }


    ## counting
    for (i in seq_len(n_cells)) {
      ii <- cells[i]
      jj <- ceiling(i / 500)

      f1_sub <- f1k[[jj]][f1k[[jj]]$cell_barcode == ii, ]

      ## count peaks
      out_summ[[ii]] <- count_peaks(
        peak_sets = peak_sets,
        filtered_fragments = f1_sub,
        extend_size = extend_size,
        n_features = n_features
      )
      # progress
      if (verbose) {
        pb$tick()
      }
    }

    if (verbose) {
      cat("Summarizing cell-by-peak matrix\n")
    }

    ## convert to a sparse matrix
    out_mat <- list_to_sparseMatrix(
      list_s_vetors = out_summ,
      n_features = n_features
    )

    ## add peak information into rownames of output
    rownames(out_mat) <- .format_peak_names(peak_sets)
  } else {
    ## use Rsamtools to load data
    tbx <- Rsamtools::TabixFile(fragment_tsv_gz_file_location)
    ## print job status
    if (verbose) {
      cat("Data loaded by chromosome\n")
    }

    ## Group peaks by chromosome while retaining their original indices. This
    ## lets us count chromosome-local matrices and restore arbitrary input order.
    peak_seqnames <- as.character(GenomeInfoDb::seqnames(peak_sets))
    slevels <- unique(peak_seqnames)
    peak_indices <- split(
      seq_len(n_features),
      factor(peak_seqnames, levels = slevels)
    )
    tabix_seqnames <- Rsamtools::seqnamesTabix(tbx)

    ## save final output
    out_mat_seq <- vector("list", length(slevels))
    names(out_mat_seq) <- slevels

    if (verbose) {
      cat("Computing peak vector for each cell.\n")
      pb <- progress::progress_bar$new(
        total = length(slevels),
        format = "[:bar] :percent computed, elapsed: :elapsed",
        clear = FALSE,
        width = 60
      )
    }

    ## load data for each chromosome
    for (sind in seq_along(slevels)) {
      seq_name <- slevels[sind]
      seq_peak_indices <- peak_indices[[seq_name]]
      seq_peak_sets <- peak_sets[seq_peak_indices]
      n_features_seq <- length(seq_peak_sets)

      query_padding <- ceiling(extend_size / 2)
      query <- GenomicRanges::GRanges(
        seqnames = seq_name,
        ranges = IRanges::IRanges(
          start = max(1L, min(BiocGenerics::start(seq_peak_sets)) - query_padding),
          end = max(BiocGenerics::end(seq_peak_sets)) + query_padding
        )
      )

      tabix_lines <- character()
      if (seq_name %in% tabix_seqnames) {
        tabix_lines <- Rsamtools::scanTabix(tbx, param = query)[[1L]]
      }

      if (length(tabix_lines) == 0L) {
        f1_seq <- GenomicRanges::GRanges()
      } else {
        f1_seq <- data.table::fread(
          text = paste(tabix_lines, collapse = "\n"),
          sep = "\t",
          header = FALSE,
          select = 1:4,
          showProgress = FALSE
        )
        f1_seq <- as.data.frame(f1_seq)
        colnames(f1_seq) <- c("seqname", "start", "end", "cell_barcode")
        f1_seq <- f1_seq[f1_seq$cell_barcode %in% cells, , drop = FALSE]
        f1_seq <- GenomicRanges::makeGRangesFromDataFrame(
          f1_seq,
          keep.extra.columns = TRUE
        )
      }

      if (verbose) {
        pb$tick()
      }

      ## create temporal output object for each seqlevels
      out_summ <- vector("list", length = n_cells)
      names(out_summ) <- cells

      zero_vec <- Matrix::sparseVector(
        i = integer(),
        x = numeric(),
        length = n_features_seq
      )
      ## count for each cell
      for (i in seq_len(n_cells)) {
        ii <- cells[i]
        f1_sub <- f1_seq[f1_seq$cell_barcode == ii, ]

        if (length(f1_sub) == 0) {
          out_summ[[ii]] <- zero_vec
          next
        }

        ## deduplicate f1_sub
        if (deduplicate) {
          f1_sub <- unique(f1_sub)
        }

        out_summ[[ii]] <- count_peaks(
          peak_sets = seq_peak_sets,
          filtered_fragments = f1_sub,
          extend_size = extend_size,
          n_features = n_features_seq
        )
      }
      if (verbose) {
        cat("Summarizing cell-by-peak matrix\n")
      }

      ## convert to a sparse matrix
      out_mat_seq[[seq_name]] <- list_to_sparseMatrix(
        list_s_vetors = out_summ,
        n_features = n_features_seq
      )

      rownames(out_mat_seq[[seq_name]]) <- .format_peak_names(seq_peak_sets)
    }

    out_mat <- do.call(rbind, out_mat_seq)
    assembled_indices <- unlist(peak_indices, use.names = FALSE)
    out_mat <- out_mat[order(assembled_indices), , drop = FALSE]
    rownames(out_mat) <- .format_peak_names(peak_sets)
  }
  out_mat
}
