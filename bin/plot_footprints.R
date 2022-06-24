

.groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gs <- lapply(unique(groups), function(x){
    if (sparse){
      matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
    }else{
      matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gs) <- unique(groups)
  return(gs)
}

.groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

.ggFootprint <- function(
  seFoot = NULL,
  name = NULL,
  pal = NULL,
  smoothWindow = NULL,
  flank = NULL,
  flankNorm = NULL,
  baseSize = 6,
  normMethod = NULL,
  logFile = NULL
  ){

  errorList <- list()

  #Get Footprint Info
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="footprint"),], name)
  biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="bias"),], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
  biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

  errorList$footMat <- footMat
  errorList$biasMat <- biasMat
  errorList$footDF <- footDF
  errorList$biasDF <- biasDF

  #Smooth Foot and Bias Mat because of sparsity
  if(!is.null(smoothWindow)){
    # .logMessage("Applying smoothing window to footprint", logFile = logFile)
    footMat <- apply(footMat, 2, function(x) .centerRollMean(x, smoothWindow))
    biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x, smoothWindow))
  }

  #Normalize Foot and Bias Mat
  # .logMessage("Normalizing by flanking regions", logFile = logFile)
  idx <- which(abs(footDF$x) >= flank - flankNorm)
  footMat <- t(t(footMat) / colMeans(footMat[idx, ,drop=FALSE]))
  biasMat <- t(t(biasMat) / colMeans(biasMat[idx, ,drop=FALSE]))

  errorList$footMatNorm <- footMat
  errorList$biasMatNorm <- footMat

  #Norm Foot By Bias
  if(tolower(normMethod) == "none"){
    title <- ""
  }else if(tolower(normMethod) == "subtract"){
    title <- "Tn5 Bias Subtracted\n"
    footMat <- footMat - biasMat
  }else if(tolower(normMethod) == "divide"){
    title <- "Tn5 Bias Divided\n"
    footMat <- footMat / biasMat
  }else{
    stop("normMethod not recognized!")
  }
  # .logMessage(paste0("NormMethod = ", normMethod), logFile = logFile)

  #Get Mean and SD for each Assay
  footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
  footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) .centerRollMean(x, 11)))

  errorList$footMatMean <- footMatMean
  errorList$footMatSd <- footMatSd
  errorList$biasMatMean <- biasMatMean
  errorList$biasMatSd <- biasMatSd
  errorList$smoothFoot <- smoothFoot

  #Create Plot Data Frames
  plotIdx <- seq_len(nrow(footMatMean)) #sort(unique(c(1, seq(1, nrow(footMatMean), smoothWindow), nrow(footMatMean))))
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x){
    data.frame(
      x = footDF$x,
      mean = footMatMean[,x],
      sd = footMatSd[,x],
      group = colnames(footMatMean)[x]
      )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))

  plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x){
    data.frame(
      x = biasDF$x,
      mean = biasMatMean[,x],
      sd = biasMatSd[,x],
      group = colnames(biasMatMean)[x]
      )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))

  errorList$plotFootDF <- plotFootDF
  errorList$plotBiasDF <- plotBiasDF

  out <- tryCatch({

    #Plot GG
    if(is.null(pal)){
      pal <- paletteDiscrete(values=gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
    }

    plotMax <- plotFootDF[order(plotFootDF$mean,decreasing=TRUE),]
    plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 50, ] #<= flank - flankNorm,]
    plotMax <- plotMax[!duplicated(plotMax$group),]
    plotMax <- plotMax[seq_len(ceiling(nrow(plotMax) / 4)), ]
    plotMax$x <- 25

    ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) +
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
      geom_line() +
      ggrepel::geom_text_repel(data = plotMax, aes(label = group), size = 3, xlim = c(75, NA)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      xlab("Distance to motif center (bp)") +
      coord_cartesian(
        expand = FALSE,
        ylim = c(quantile(plotFootDF$mean, 0.0001), 1.15*quantile(smoothFoot, 0.999)),
        xlim = c(min(plotFootDF$x),max(plotFootDF$x))
      ) + theme_ArchR(baseSize = baseSize) + ggtitle(name) +
      guides(fill = FALSE) +
      guides(color = FALSE) + ylab(paste0(title,"Normalized Insertions"))
      # ggrepel::geom_label_repel(data = plotMax, aes(label = group), size = 3, xlim = c(75, NA))

    ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) +
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
      geom_line() +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      xlab("Distance to motif center (bp)") +
      coord_cartesian(
        expand = FALSE,
        ylim = c(quantile(plotBiasDF$mean, 0.0001), 1.05*quantile(plotBiasDF$mean, 0.999)),
        xlim = c(min(plotBiasDF$x),max(plotBiasDF$x))
      ) + theme_ArchR(baseSize = baseSize) + ylab("Tn5-Bias Normalized Insertions") +
      theme(legend.position = "bottom", legend.box.background = element_rect(color = NA))

    ggAlignPlots(ggFoot, .ggSmallLegend(ggBias), sizes=c(2,1), draw = FALSE)

  }, error = function(e){

    #.logError(e, fn = ".ggFootprint", info = name, errorList = errorList, logFile = logFile)
    print(e)
  })

  out

}

.ggSmallLegend <- function(
  gg = NULL,
  pointSize = 2,
  baseSize = 5,
  spaceLegend = 0.1
  ) {
    #https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
    gg +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = baseSize),
              legend.text  = element_text(size = baseSize),
              legend.key.size = unit(spaceLegend, "lines"))
}

.centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}


.getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}

plotFootprintsManual <- function(
  seFoot = NULL,
  names = NULL,
  pal = NULL,
  flank = 250,
  flankNorm = 50,
  normMethod = "Subtract",
  smoothWindow = NULL,
  baseSize = 6,
  plot = TRUE,
  ArchRProj = NULL,
  plotName = paste0("Plot-Footprints-", normMethod),
  height = 6,
  width = 4,
  addDOC = TRUE,
  force = FALSE,
  logFile = createLogFile("plotFootprints")
  ){

  tstart <- Sys.time()

  if(is.null(names)){
    names <- names(assays(seFoot))
  }


  if(plot){

    name <- gsub("\\.pdf", "", plotName)
    if(is.null(ArchRProj)){
      outDir <- "Plots"
    }else{
      # ArchRProj <- .validArchRProject(ArchRProj)
      outDir <- file.path(getOutputDirectory(ArchRProj), "Plots")
    }

    dir.create(outDir, showWarnings = FALSE)
    if(addDOC){
      doc <- gsub(":","-",stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2])
      filename <- file.path(outDir, paste0(name, "_Date-", Sys.Date(), "_Time-", doc, ".pdf"))
    }else{
      filename <- file.path(outDir, paste0(name, ".pdf"))
    }

    pdf(filename, width = width, height = height, useDingbats = FALSE)

  }

  ggList <- lapply(seq_along(names), function(x){


    gg <- .ggFootprint(
      seFoot = seFoot,
      name = names[x],
      pal = pal,
      smoothWindow = smoothWindow,
      flank = flank,
      flankNorm = flankNorm,
      baseSize = baseSize,
      normMethod = normMethod,
      logFile = logFile
    )

    if(plot){
      if(x != 1){
        grid::grid.newpage()
      }
      grid::grid.draw(gg)
      return(0)
    }else{
      return(gg)
    }

  })
  

  if(!plot){
    names(ggList) <- names
    ggList
  }else{
    dev.off()
    return(invisible(0))
  }
}
