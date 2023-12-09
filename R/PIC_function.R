## PIC function



#' Function for loading fragment files and filter by cell barcodes
#'
#' @param fragment_tsv_gz_file_location The 10X Cell Ranger output
#'  fragment.tsv.gz file location. This can usually be found at the /out
#'  directory from Cell Ranger output
#' @param cells The cell barcode lables as a Character vector
#' @param verbose Whether to output progress message. Default TRUE
#'
#' @return data.frame containing fragments that are filtered by cell barcodes
#' @export
load_fragments <- function(
    fragment_tsv_gz_file_location, cells, verbose = TRUE) {
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
  n_cells_fragment_file <- sum(cells_retain)

  ## report the proportion of reads in cell barcodes
  prop_bcreads <- n_cells_fragment_file / dim(f1)[1]
  if (verbose) {
    print(paste("proportion of reads in cell barcdes is ", prop_bcreads,
      sep = ""
    ))
  }

  ## error when no cells found in the fragment file
  if (n_cells_fragment_file < 1) {
    stop("Cell barcodes not found in fragment files, please check input")
  } else if (n_cells_fragment_file < 10) {
    warning("Fewer than 10 cells found in fragment files,
                please consider checking input")
  }

  f1 <- f1[cells_retain, ]
  return(f1)
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
                        extend_size, n_features) {
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

  out_vec <- sparseVector(x = counts$Freq, i = counts$ol, length = n_features)

  return(out_vec)
}


#' Generate a sparse matrix from a list of sparse vectors
#'
#' @param chunk A chunk of lists with each element being a sparseVector
#'
#' @return A sparseMatrix that is a column bind of all sparseVectors within
#'  the chunk
#' @noRd
make_s_mat_from_s_vec <- function(chunk) {
  indices <- lapply(chunk, function(y) y@i)
  values <- lapply(chunk, function(y) y@x)

  i <- unlist(indices)
  j <- unlist(lapply(
    seq_along(chunk),
    function(k) rep(k, length(indices[[k]]))
  ))
  x <- unlist(values)

  s_mat <- sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(chunk[[1]]), length(chunk))
  )

  ## add cell names to the matrix
  colnames(s_mat) <- names(chunk)
  return(s_mat)
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
  n_cells <- length(list_s_vetors)
  chunk_size <- ceiling(n_cells * n_features / 2^31)

  if (chunk_size > 1) {
    chunks <- split(
      list_s_vetors,
      ceiling(seq_along(list_s_vetors) / chunk_size)
    )
    sparse_matrices <- lapply(chunks, make_s_mat_from_s_vec)
    sparse_matrices <- do.call(cbind, sparse_matrices)
  } else {
    sparse_matrices <- make_s_mat_from_s_vec(list_s_vetors)
  }
  return(sparse_matrices)
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
  ## if colnames not specified, we specify by order
  if (is.null(colnames(peak_sets)) ||
    !(all(c("seqname", "start", "end") %in% colnames(peak_sets)))) {
    colnames(peak_sets) <- c("seqname", "start", "end")
  }

  ## convert into GRanges
  peak_sets <- try(GenomicRanges::makeGRangesFromDataFrame(peak_sets),
    silent = TRUE
  )
  if (inherits(peak_sets, "try-error")) {
    stop("An error occurred in trying to convert peak_sets
          into a GRanges object, please check input")
  } else {
    return(peak_sets)
  }
}



#' Title PIC-counting data matrix
#'
#' @importFrom methods is
#' @importFrom data.table fread
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @import Rsamtools
#' @importFrom IRanges subsetByOverlaps
#' @import Matrix
#' @importFrom utils read.csv str
#' @import progress
#' @param cells The cell barcode lables as a Character vector
#' @param fragment_tsv_gz_file_location The 10X Cell Ranger output
#'  fragment.tsv.gz file location. This can usually be found at the /out
#'  directory from Cell Ranger output
#' @param peak_sets The set of peaks as a GenomicRanges object. This will be
#'  the features for the data matrix. Alternatively, this can be a data.frame
#'  and the function will convert it into a GenomicRanges object
#' @param deduplicate Whether to include deduplicate step where within
#'  the same cell,
#'  fragments with identical start and end location will be deduplicated.
#'  This is usually unnecessisary from Cell Ranger ATAC output, since
#'  Cell Ranger ATAC has already deduplicated the fragments.
#'  But for dsc-ATAC-seq data, this step will
#'  be helpful and recommended.
#' @param load_full Whether to load the whole fragment.tsv.gz file into memory.
#'  If set to FALSE, the function will load it dynamically to save RAM
#' @param extend_size How long should we extend the exact insertion site as
#'  accessible window
#' @param verbose Whether to output progress information including the progress
#'  bar
#'
#' @return The peak by cell PIC count matrix
#' @export
#'
PIC_counting <- function(cells,
                         fragment_tsv_gz_file_location,
                         peak_sets,
                         deduplicate = FALSE,
                         load_full = TRUE,
                         extend_size = 5L,
                         verbose = TRUE) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The GenomicRanges package is not installed.
         Please install it using BiocManager::install('GenomicRanges').")
  }

  ## check input
  if (extend_size < 0) {
    stop("extend_size has to be a positive integer!")
  }

  n_cells <- length(cells)
  n_features <- length(peak_sets)

  if (n_cells == 0 | anyNA(cells)) {
    stop("cell names are empty or contain NA values!")
  }

  ## we accpt peak_sets to be a GRanges or we convert it into one
  if (!methods::is(peak_sets, "GRanges")) {
    peak_sets <- data_frame_to_GRanges(peak_sets)
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
    if (sum(f1$start - 1 >= f1$end) >= 1) {
      f1_sub1 <- f1[f1$start - 1 < f1$end, ]
      f1_sub2 <- f1[f1$start - 1 >= f1$end, ]
      f1_sub2s <- f1_sub2$start
      f1_sub2$start <- f1_sub2$end
      f1_sub2$end <- f1_sub2s
      f1 <- rbind(f1_sub1, f1_sub2)
      rm(f1_sub2)
      rm(f1_sub2s)
      rm(f1_sub1)
    }


    ## generate GenomicRanges object
    f1 <- GenomicRanges::makeGRangesFromDataFrame(f1,
      keep.extra.columns = T
    )

    ## pre-sort fragments
    f1_s <- IRanges::subsetByOverlaps(f1,
      ranges = peak_sets,
      maxgap = ceiling(extend_size / 2)
    )
    rm(f1)

    n_subset <- n_cells %/% 500 + 1
    f1k <- rep(list(), length = n_subset)
    for (i in 1:n_subset) {
      s <- (i - 1) * 500 + 1
      e <- min(i * 500, n_cells)

      ## deduplicate f1_s
      if (deduplicate) {
        f1k[[i]] <- unique(f1_s[f1_s$cell_barcode %in% cells[s:e], ])
      }else{
        f1k[[i]] <- f1_s[f1_s$cell_barcode %in% cells[s:e], ]
      }

    }
    rm(f1_s)
    gc()

    ## progress bar
    pb <- progress::progress_bar$new(
      total = n_cells,
      format = "computing peak vector for each cell",
      clear = FALSE
    )

    ## counting
    for (i in 1:n_cells) {
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
      pb$tick()
    }

    ## convert to a sparse matrix
    out_mat <- list_to_sparseMatrix(
      list_s_vetors = out_summ,
      n_features = n_features
    )

    ## add peak information into rownames of output
    pkdf <- as.data.frame(peak_sets)
    fname <- paste(pkdf$seqnames, ":", pkdf$start, "-", pkdf$end, sep = "")
    rownames(out_mat) <- fname
  } else {
    ## use Rsamtools to load data
    tbx <- Rsamtools::TabixFile(fragment_tsv_gz_file_location)

    ## get ranges for each chromosome
    slevels <- seqlevels(peak_sets)
    grl <- split(peak_sets, seqnames(peak_sets))
    grl <- sort(grl)
    n_features_seq <- sapply(grl, length)
    param <- matrix(nrow = length(slevels), ncol = 3)
    colnames(param) <- c("seqname", "start", "end")
    rownames(param) <- slevels
    param <- as.data.frame(param)
    param$seqname <- slevels

    for (sl in slevels) {
      param[sl, "start"] <- min(start(grl[[sl]]))
      param[sl, "end"] <- max(end(grl[[sl]]))
    }
    param <- GenomicRanges::makeGRangesFromDataFrame(param)

    ## save final output
    out_mat_seq <- rep(list(), length = length(slevels))
    names(out_mat_seq) <- slevels

    ## print job status
    print("loading data for each chromosome")

    ## progress bar
    pb <- progress::progress_bar$new(
      total = length(slevels),
      format = "computing peak vector for each chromosome",
      clear = FALSE
    )

    ## load data for each chromosome
    for (sind in seq_along(slevels)) {
      seq_name <- slevels[sind]
      res <- Rsamtools::scanTabix(tbx, param = param[sind])
      # length(res[[1]])
      f1_seq <- read.csv(textConnection(res[[1]]), sep = "\t", header = FALSE)
      f1_seq <- f1_seq[, 1:4]
      colnames(f1_seq) <- c("seqname", "start", "end", "cell_barcode")
      f1_seq <- f1_seq[f1_seq$cell_barcode %in% cells, ]
      f1_seq <- GenomicRanges::makeGRangesFromDataFrame(f1_seq,
        keep.extra.columns = TRUE
      )

      pb$tick()

      ## create temporal output object for each seqlevels
      out_summ <- rep(list(), length = n_cells)
      names(out_summ) <- cells

      zero_vec <- rep(0, length = n_features_seq[seq_name])
      ## count for each cell
      for (i in 1:n_cells) {
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
          peak_sets = peak_sets,
          filtered_fragments = f1_sub,
          extend_size = extend_size,
          n_features = n_features
        )
      }

      ## convert to a sparse matrix
      out_mat_seq[[seq_name]] <- list_to_sparseMatrix(
        list_s_vetors = out_summ,
        n_features = n_features
      )

      pkdf <- as.data.frame(grl[[seq_name]])
      fname <- paste(pkdf$seqnames, ":", pkdf$start, "-", pkdf$end, sep = "")
      rownames(out_mat_seq[[seq_name]]) <- fname
    }

    out_mat <- do.call(rbind, out_mat_seq)
    pkdf_full <- as.data.frame(peak_sets)
    fname_full <- paste(pkdf_full$seqnames, ":", pkdf_full$start,
      "-", pkdf_full$end,
      sep = ""
    )
    outmat <- outmat[fname_full, ]
  }
  return(out_mat)
}
