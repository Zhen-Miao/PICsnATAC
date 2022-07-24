## PIC_function in ArchR


#' Add a Peak Matrix to the ArrowFiles of an ArchRProject
#' 
#' This function, for each sample, will independently compute counts for each peak
#' per cell in the provided ArchRProject using the "PeakMatrix".
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param ceiling The maximum counts per feature allowed. This is used to prevent large biases in peak counts.
#' @param binarize A boolean value indicating whether the peak matrix should be binarized prior to storage. This can be useful
#' for downstream analyses when working with insertion counts.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param threads The number of threads to be used for parallel computing.
#' @param parallelParam A list of parameters to be passed for biocparallel/batchtools parallel computing.
#' @param force A boolean value indicating whether to force the "PeakMatrix" to be overwritten if it already exist in the given `ArchRProject`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addPeakMatrix_PIC <- function(
  ArchRProj = NULL,
  PIC_mat,
  ceiling = 4, 
  binarize = FALSE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = ceiling, name = "ceiling", valid = c("numeric"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  if(is.null(ArchRProj@peakSet)){
    stop("No peakSet found in ArchRProject!")
  }
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  allCells <- rownames(getCellColData(ArchRProj))
  outDir <- getOutputDirectory(ArchRProj)
  
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }
  
  .startLogging(logFile = logFile)
  .logThis(ArchRProj@peakSet, "peakSet", logFile = logFile)
  
  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$matrixName = "PeakMatrix"
  args$features <- ArchRProj@peakSet
  args$X <- seq_along(ArrowFiles)
  args$FUN <- .addFeatureMatrix
  args$registryDir <- file.path(outDir, "CountPeaksRegistry")
  
  #Remove project from args
  args$ArchRProj <- NULL
  
  #Run With Parallel or lapply
  outList <- .batchlapply(args)
  
  readsInPeaks <- lapply(outList, function(x) x$RIP) %>% unlist
  FRIP <- lapply(outList, function(x) x$FRIP) %>% unlist
  ArchRProj <- addCellColData(ArchRProj, data = readsInPeaks, name = "ReadsInPeaks", names(readsInPeaks), force = force)
  ArchRProj <- addCellColData(ArchRProj, data = FRIP, name = "FRIP", names(readsInPeaks), force = force)
  
  .endLogging(logFile = logFile)
  
  return(ArchRProj)
  
}

