## PIC function


#' Title PIC-count
#'
#' @importFrom rtracklayer import
#' @importFrom methods as
#' @import dplyr
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps
#' @import Matrix
#' @param cells The cell barcode lables as a Character vector
#' @param fragment_tsv_gz_file_location The 10X Cell Ranger output fragment.tsv.gz file location
#' @param peak_sets The set of peaks as a GenomicRanges object
#'
#' @return The counted cell by peak matrix
#' @export
#'
#' @examples
PIC_counting <- function(cells,
                         fragment_tsv_gz_file_location,
                         peak_sets) {

  ## load fragment files
  f1 <- data.table::fread(fragment_tsv_gz_file_location,header = F, 
                          select = 1:4)
  colnames(f1) <- c('seqname','start','end','name')
  f1 <- f1[name %in% cells]
  f1 = GenomicRanges::makeGRangesFromDataFrame(f1)
  
  ## create the object to save output
  n_cells <- length(cells)
  n_features <- length(peak_sets)
  out_summ <- rep(list(), length = n_cells)
  names(out_summ) <- cells

  ## pre-sort fragments
  f1_s <- subsetByOverlaps(f1, ranges = peak_sets)
  rm(f1)
  gc()
  n_subset <- n_cells %/% 500 + 1
  f1k <- rep(list(), length = n_subset)
  for (i in 1:n_subset) {
    s <- (i - 1) * 500 + 1
    e <- min(i * 500, n_cells)
    f1k[[i]] <- f1_s[f1_s$name %in% cells[s:e], ]
  }
  rm(f1_s)
  gc()


  ## counting
  for (i in 1:n_cells) {
    ii <- cells[i]
    if (i %% 500 == 0) {
      jj <- i %/% 500
    } else {
      jj <- (i %/% 500) + 1
    }

    f1_sub <- f1k[[jj]][f1k[[jj]]$name == ii, ]
    
    ## deduplicate f1_sub
    f1_sub <- unique(f1_sub)

    ## get start and end position
    f1s <- GenomicRanges::resize(f1_sub, width = 1, fix = "start")
    f1e <- GenomicRanges::resize(f1_sub, width = 1, fix = "end")

    overlaped_s <- GenomicRanges::findOverlaps(f1s, peak_sets, select = "first")
    overlaped_e <- GenomicRanges::findOverlaps(f1e, peak_sets, select = "first")
    overlaped_e_nodc <- overlaped_e

    overlaped_e_nodc[overlaped_s == overlaped_e_nodc] <- NA
    ol <- c(overlaped_s, overlaped_e_nodc)

    counts <- as.data.frame(table(ol), stringsAsFactors = F)
    counts$ol <- as.integer(counts$ol)

    out_summ[[ii]] <- sparseVector(x = counts$Freq, i = counts$ol, length = n_features)
    # print(i)
  }

  ## convert to a sparse matrix
  out_summ2 <- lapply(out_summ, as, "vector")

  n_divide <- ceiling(n_cells / 4)
  out_mat11 <- do.call(cbind, out_summ2[1:n_divide])
  out_mat12 <- do.call(cbind, out_summ2[(n_divide + 1):(n_divide * 2)])
  out_mat13 <- do.call(cbind, out_summ2[(n_divide * 2 + 1):(n_divide * 3)])
  out_mat14 <- do.call(cbind, out_summ2[(n_divide * 3 + 1):n_cells])

  o1 <- as(out_mat11, "sparseMatrix")
  o2 <- as(out_mat12, "sparseMatrix")
  o3 <- as(out_mat13, "sparseMatrix")
  o4 <- as(out_mat14, "sparseMatrix")

  out_mat <- cbind(o1, o2, o3, o4)
  pkdf <- as.data.frame(peak_sets)
  fname <- paste(pkdf$seqnames, ":", pkdf$start, "-", pkdf$end, sep = "")
  rownames(out_mat) <- fname

  return(out_mat)
}
