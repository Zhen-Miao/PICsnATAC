## PIC function


#' Title PIC-count
#'
#' @importFrom rtracklayer import
#' @importFrom methods as
#' @import dplyr
#' @import data.table
#' @import GenomicRanges
#' @import Rsamtools
#' @importFrom IRanges subsetByOverlaps
#' @import Matrix
#' @param cells The cell barcode lables as a Character vector
#' @param fragment_tsv_gz_file_location The 10X Cell Ranger output fragment.tsv.gz file location
#' @param peak_sets The set of peaks as a GenomicRanges object
#' @param deduplicate Whether to include deduplicate step where within the same cell,
#'  fragments with identical start and end location will be deduplicated. This is
#'  usually unnecessisary from Cell Ranger ATAC output, since Cell Ranger ATAC has
#'  already deduplicated the fragments. But for dsc-ATAC-seq data, this step will 
#'  be helpful and recommended.
#' @param load_full Whether to load the whole fragment.tsv.gz file into memory. If set
#'  to FALSE, the function will load it dynamically to save RAM
#'
#' @return The counted cell by peak matrix
#' @export
#'
#' @examples
PIC_counting <- function(cells,
                         fragment_tsv_gz_file_location,
                         peak_sets,
                         deduplicate = FALSE,
                         load_full = TRUE) {
  
  ## create the object to save output
  n_cells <- length(cells)
  n_features <- length(peak_sets)
  
  ## if load full files
  if(load_full){
    ## create output object
    out_summ <- rep(list(), length = n_cells)
    names(out_summ) <- cells
    
    ## load fragment files
    f1 <- data.table::fread(fragment_tsv_gz_file_location,header = F, 
                            select = 1:4)
    # colnames(f1) <- c('seqname','start','end','cell_barcode')
    setnames(f1, c('seqname','start','end','cell_barcode'))
    f1 <- f1[cell_barcode %in% cells]
    f1 = GenomicRanges::makeGRangesFromDataFrame(f1,
                                                 keep.extra.columns=T)
    
    ## pre-sort fragments
    f1_s <- subsetByOverlaps(f1, ranges = peak_sets)
    rm(f1)
    gc()
    n_subset <- n_cells %/% 500 + 1
    f1k <- rep(list(), length = n_subset)
    for (i in 1:n_subset) {
      s <- (i - 1) * 500 + 1
      e <- min(i * 500, n_cells)
      f1k[[i]] <- f1_s[f1_s$cell_barcode %in% cells[s:e], ]
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
      
      f1_sub <- f1k[[jj]][f1k[[jj]]$cell_barcode == ii, ]
      
      ## deduplicate f1_sub
      if(deduplicate){
        f1_sub <- unique(f1_sub)
      }
      
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
    

  } else {
    ## use Rsamtools to load data
    tbx <- TabixFile(fl)
    
    ## get ranges for each chromosome
    slevels <- seqlevels(peak_sets)
    grl <- split(peak_sets, seqnames(peak_sets))
    grl <- sort(grl)
    n_features_seq <- sapply(grl, length)
    param <- matrix(nrow = length(slevels),ncol = 3)
    colnames(param) <- c('seqname', 'start', 'end')
    rownames(param) <- slevels
    param <- as.data.frame(param)
    param$seqname <- slevels
    
    for(sl in slevels){
      param[sl,'start'] <- min(start(grl[[sl]]))
      param[sl,'end'] <- max(end(grl[[sl]]))
    }
    param <- GenomicRanges::makeGRangesFromDataFrame(param)
    
    ## save final output
    out_mat_seq <- rep(list(), length = length(slevels))
    names(out_mat_seq) <- slevels
    
    ## load data for each chromosome
    for(sind in seq_along(slevels)){ 
      seq_name <- slevels[sind]
      res <- scanTabix(tbx, param=param[sind])
      # length(res[[1]])
      f1_seq <- read.csv(textConnection(res[[1]]), sep="\t", header=FALSE)
      f1_seq <- f1_seq[,1:4]
      colnames(f1_seq) <- c('seqname','start','end','cell_barcode')
      f1_seq <- f1_seq[f1_seq$cell_barcode %in% cells,]
      f1_seq = GenomicRanges::makeGRangesFromDataFrame(f1_seq,
                                                       keep.extra.columns=T)
      
      ## create temporal output object for each seqlevels 
      out_summ <- rep(list(), length = n_cells)
      names(out_summ) <- cells
      
      zero_vec <- rep(0,length = n_features_seq[seq_name])
      ## count for each cell
      for (i in 1:n_cells) {
        
        ii <- cells[i]
        f1_sub <- f1_seq[f1_seq$cell_barcode == ii, ]
        
        if(length(f1_sub) == 0){
          out_summ[[ii]] <- zero_vec
          next
        }
        
        ## deduplicate f1_sub
        if(deduplicate){
          f1_sub <- unique(f1_sub)
        }
        
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
        
        out_summ[[ii]] <- sparseVector(x = counts$Freq, i = counts$ol, 
                                       length = n_features_seq[seq_name])
        # print(i)
      }
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
      
      out_mat_seq[[seq_name]] <- cbind(o1, o2, o3, o4)
      pkdf <- as.data.frame(grl[[seq_name]])
      fname <- paste(pkdf$seqnames, ":", pkdf$start, "-", pkdf$end, sep = "")
      rownames(out_mat_seq[[seq_name]]) <- fname
    }
    out_mat <- do.call(rbind, out_mat_seq)
    pkdf_full <- as.data.frame(peak_sets)
    fname_full <- paste(pkdf_full$seqnames, ":", pkdf_full$start, "-", pkdf_full$end,
                   sep = "")
    outmat <- outmat[fname_full,]
    
  }
  return(out_mat)
}
