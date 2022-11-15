## all other methods for differential test
library(Rfast)
library(presto)
library(SummarizedExperiment)

## the same thing, but accelerated 
seurat_method2_acc <- function(data_matrix_pos, 
                               data_matrix_neg,
                               peak_region_fragments){
  data.use <- cbind(data_matrix_pos, data_matrix_neg)
  group.info <- c(rep(0,times = dim(data_matrix_pos)[2]),
                  rep(1,times = dim(data_matrix_neg)[2]) )
  
  # group.info <- as.factor(group.info)
  # peak_region_fragments <- Rfast::colsums(data.use)
  ## data.use should be feature by cell matrix
  n_features <- dim(data.use)[1]
  p_val <- vector(length = n_features)
  for(i in 1:n_features){
    X1 = cbind(data.use[i, ],peak_region_fragments)
    # model.data <- data.frame(GENE = data.use[i, ], group.info = group.info,
    #                          peak_region_fragments = peak_region_fragments)
    # model1 <- glm(formula = group.info ~ GENE + peak_region_fragments,
    #               data = model.data, family = binomial)
    # model2 <- glm(formula = group.info ~ peak_region_fragments,
    #               data = model.data, family = binomial)
    model1 <- glm_logistic(x = X1, y = group.info)
    model2 <- glm_logistic(x = peak_region_fragments,y = group.info)
    p_val[i] <- pchisq(model2$devi - model1$devi, df = 1,lower.tail = F)
  }
  return(p_val)
}


## archR functionality
"%ni%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)


getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}

computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
){
  # .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  # .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  # .validInput(input = k, name = "k", valid = c("integer"))
  # .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  # .requirePackage("nabor", source = "cran")
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


matchBiasCellGroups <- function(
  input = NULL,
  groups = NULL,
  useGroups = NULL,
  bgdGroups = NULL,
  bias = NULL,
  k = 100,
  n = 500,
  bufferRatio = 0.8
){
  
  #Summary Function
  .summarizeColStats <- function(m = NULL, name = NULL){
    med <- apply(m, 2, median)
    mean <- colMeans(m)
    sd <- apply(m, 2, sd)
    loQ <- apply(m, 2, function(x) quantile(x, 0.25))
    hiQ <- apply(m, 2, function(x) quantile(x, 0.75))
    summaryDF <- t(data.frame(
      median = med,
      mean = mean,
      sd = sd,
      lowerQuartile = loQ,
      upperQuartile = hiQ
    )) %>% data.frame
    colnames(summaryDF) <- colnames(m)
    if(!is.null(name)){
      summaryDF$name <- name
    }
    summaryDF
  }
  

  
  #Make sure input is dataframe
  input <- data.frame(input)
  
  #Norm using input string ie log10(nfrags)
  inputNorm <- lapply(seq_along(bias), function(x){
    plyr::mutate(input, o=eval(parse(text=bias[x])))$o
  }) %>% Reduce("cbind", .) %>% data.frame
  rownames(inputNorm) <- rownames(input)
  
  #Quantile Normalization
  inputNormQ <- lapply(seq_len(ncol(inputNorm)), function(x){
    getQuantiles(inputNorm[,x])
  }) %>% Reduce("cbind", .) %>% data.frame
  rownames(inputNormQ) <- rownames(input)
  
  #Add Colnames
  colnames(inputNorm) <- bias
  colnames(inputNormQ) <- bias
  
  if(is.null(useGroups)){
    useGroups <- gtools::mixedsort(unique(paste0(groups)))
  }
  
  if(is.null(bgdGroups)){
    bgdGroups <- gtools::mixedsort(unique(paste0(groups)))
  }
  
  stopifnot(all(useGroups %in% unique(paste0(groups))))
  stopifnot(all(bgdGroups %in% unique(paste0(groups))))
  
  #Get proportion of each group
  prob <- table(groups) / length(groups)
  bgdProb <- prob[which(names(prob) %in% bgdGroups)] / sum(prob[which(names(prob) %in% bgdGroups)])
  
  #pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  matchList <- lapply(seq_along(useGroups), function(x){
    
    #setTxtProgressBar(pb,round(x*100/length(useGroups),0))
    
    #############
    # Organize
    #############
    groupx <- useGroups[x]
    idx <- which(names(bgdProb) == groupx)
    if(length(idx) > 0 & length(idx) != length(bgdProb)){
      bgdProbx <- bgdProb[-idx]/sum(bgdProb[-idx])
    }else{
      bgdProbx <- bgdProb
    }
    
    idF <- which(groups == groupx)
    
    if(all(length(idF) * bgdProbx < 1)){
      if(length(idF) < length(bgdProbx)){
        bgdProbx <- bgdProbx[sample(names(bgdProbx), floor(length(idF) * bufferRatio))]
        bgdProbx[1:length(bgdProbx)] <- rep(1/length(bgdProbx), length(bgdProbx))
      }
    }
    
    idB <- which(groups %in% names(bgdProbx))
    
    if(k > length(idB)){
      k2 <- length(idB)
    }else{
      k2 <- k
    }
    
    knnx <- computeKNN(inputNormQ[idB, ,drop=FALSE], inputNormQ[idF, ,drop=FALSE], k = k2)
    sx <- sample(seq_len(nrow(knnx)), nrow(knnx))
    
    minTotal <- min(n, length(sx) * bufferRatio)
    nx <- sort(floor(minTotal * bgdProbx))
    
    ###############
    # ID Matching
    ###############
    idX <- c()
    idY <- c()
    it <- 0
    
    if(any(nx <= 0)){
      nx[which(nx <= 0)] <- Inf
      nx <- sort(nx)
    }
    
    while(it < length(sx) & length(idX) < minTotal){
      
      it <- it + 1
      knnit <- knnx[sx[it],]
      groupit <- match(groups[idB][knnit],names(nx))
      selectUnique <- FALSE
      selectit <- 0
      oit <- order(groupit)
      
      while(!selectUnique){
        selectit <- selectit + 1
        itx <- which(oit==selectit)
        cellx <- knnit[itx]
        groupitx <- groupit[itx]
        if(is.infinite(nx[groupitx])){
          if(selectit == k2){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }else{
          if(cellx %ni% idY){
            selectUnique <- TRUE
          }
          if(selectit == k2){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }
      }
      
      if(!is.na(itx)){
        idX <- c(idX, sx[it])
        idY <- c(idY, cellx)
        nx[groupitx] <- nx[groupitx] - 1
        if(any(nx <= 0)){
          nx[which(nx <= 0)] <- Inf
          nx <- sort(nx)
        }
      }
      
      if(all(is.infinite(nx))){
        it <- length(sx)
      }
      
    }
    
    #####################
    # Convert Back to Normal Indexing
    #####################
    idX <- seq_len(nrow(inputNormQ))[idF][idX]
    idY <- seq_len(nrow(inputNormQ))[idB][idY]
    
    #####################
    # Matching Stats Groups
    #####################
    estbgd <- sort(floor(minTotal * bgdProbx))
    obsbgd <- rep(0, length(estbgd))
    names(obsbgd) <- names(estbgd)
    tabGroups <- table(groups[idY])
    obsbgd[names(tabGroups)] <- tabGroups
    estbgdP <- round(100 * estbgd / sum(estbgd),3)
    obsbgdP <- round(100 * obsbgd / sum(obsbgd),3)
    
    #####################
    # Matching Stats Bias Norm Values
    #####################
    forBias <- .summarizeColStats(inputNorm[idX,,drop=FALSE], name = "foreground")
    bgdBias <- .summarizeColStats(inputNorm[idY,,drop=FALSE], name = "background")
    
    out <- list(
      cells = idX, 
      bgd = idY, 
      summaryCells = forBias, 
      summaryBgd = bgdBias, 
      bgdGroups = rbind(estbgd, obsbgd),
      bgdGroupsProbs = rbind(estbgdP, obsbgdP),
      corbgdGroups = suppressWarnings(cor(estbgdP, obsbgdP)),
      n = length(sx), 
      p = it / length(sx),
      group = groupx,
      k = k2
    )
    
    return(out)
    
  }) %>% SimpleList
  names(matchList) <- useGroups
  
  outList <- SimpleList(
    matchbgd = matchList,
    info = SimpleList(
      cells = rownames(input),
      groups = groups,
      biasNorm = inputNorm,
      biasNormQ = inputNormQ
    )
  )
  
  return(outList)
  
}



archR_method <- function(data_matrix_pos_count, data_matrix_neg_count){
  
  n_cells_pos <- dim(data_matrix_pos_count)[2]
  n_cells_neg <- dim(data_matrix_neg_count)[2]
  colnames(data_matrix_pos_count) <- as.character(c(1:n_cells_pos))
  colnames(data_matrix_neg_count) <- as.character(c(1:n_cells_neg))
  n_fragments <- c(Matrix::colSums(data_matrix_pos_count),
                   Matrix::colSums(data_matrix_neg_count))
  

  ## first, match the background cells with ArchR functions
  ## all input identical to ArchR default parameters 
  colDat <- data.frame(logn_frag = log10(n_fragments))
  groups <- c( rep('A', times = n_cells_pos),
               rep('B', times = n_cells_neg))
  useGroups <- 'A'
  bgdGroups <- 'B'
  bias = "logn_frag"
  k = 100
  bufferRatio = 0.8
  maxCells = 5000
  
  matchObj <- matchBiasCellGroups(
    input = colDat, 
    groups = groups,
    useGroups = useGroups,
    bgdGroups = bgdGroups,
    bias = bias,
    k = k,
    n = maxCells,
    bufferRatio = bufferRatio
  )
  
  print('matchObj finished')
  
  matchx <- matchObj[[1]][[useGroups]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]
  
  
  ######## ------- show the length of matched cells ----
  print('length of matched cells')
  print(length(bgdx))
  print(head(bgdx))
  print(max(as.numeric(bgdx)))
  
  mat1 <- data_matrix_pos_count[,as.numeric(cellsx)]
  mat2 <- data_matrix_neg_count[,as.numeric(bgdx) - n_cells_pos]
  
  df <- wilcoxauc(cbind(data_matrix_pos_count,data_matrix_neg_count),
                  c(rep("Top", ncol(data_matrix_pos_count)),
                    rep("Bot", ncol(data_matrix_neg_count))))
  df <- df[which(df$group=="Top"),]
  return(df$pval)
}



archR_method_custom <- function(data_matrix_pos_count,
                                data_matrix_neg_count,
                                k = 20,
                                bufferRatio = 0.5,
                                maxCells = 5000){
  
  n_cells_pos <- dim(data_matrix_pos_count)[2]
  n_cells_neg <- dim(data_matrix_neg_count)[2]
  colnames(data_matrix_pos_count) <- as.character(c(1:n_cells_pos))
  colnames(data_matrix_neg_count) <- as.character(c(1:n_cells_neg))
  n_fragments <- c(Matrix::colSums(data_matrix_pos_count),
                   Matrix::colSums(data_matrix_neg_count))
  
  
  ## first, match the background cells with ArchR functions
  ## all input identical to ArchR default parameters 
  colDat <- data.frame(logn_frag = log10(n_fragments))
  groups <- c( rep('A', times = n_cells_pos),
               rep('B', times = n_cells_neg))
  useGroups <- 'A'
  bgdGroups <- 'B'
  bias = "logn_frag"

  
  matchObj <- matchBiasCellGroups(
    input = colDat, 
    groups = groups,
    useGroups = useGroups,
    bgdGroups = bgdGroups,
    bias = bias,
    k = k,
    n = maxCells,
    bufferRatio = bufferRatio
  )
  
  print('matchObj finished')
  
  matchx <- matchObj[[1]][[useGroups]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]
  
  
  ######## ------- show the length of matched cells ----
  print('length of matched cells')
  print(length(bgdx))
  print(head(bgdx))
  print(max(as.numeric(bgdx)))
  
  mat1 <- data_matrix_pos_count[,as.numeric(cellsx)]
  mat2 <- data_matrix_neg_count[,as.numeric(bgdx) - n_cells_pos]
  
  df <- wilcoxauc(cbind(mat1,mat2),
                  c(rep("Top", ncol(mat1)),
                    rep("Bot", ncol(mat2))))
  df <- df[which(df$group=="Top"),]
  return(df$pval)
}



library(lattice)
qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

