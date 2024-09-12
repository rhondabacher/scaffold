#' Sequencing step of the Scaffold simulation
#' @param amplifiedMolecules A list of captured molecules for each sample.
#' @param equalizationAmount A value between 0 and 1 indicating the q* to determine the number of samples that undergo dilution in the equalization step of the simulation. A value of 0 indicates all cells are diluted to the smallest concentration and a value of 1 indicates no equalization is performed.
#' @param totalDepth The total sequencing depth of the simulated data. If left NULL, this is taken from the \code{sce} object. If more cells are generated than in the original dataset, then the totalDepth will be scaled up accordingly.
#' @param efficiencyPCR A numeric vector representing the efficiency of PCR for each sample.
#' @param roundsPCR An integer indicating the number of PCR amplification rounds.
#' @param efficiencyTag A numeric vector representing the efficiency of tagging for each cell.
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param useUMI A TRUE/FALSE indicating whether the protocol should use UMIs (Unique Molecular Identifiers). Droplet or 10X protocols have this set as TRUE for the default, otherwise FALSE.
#' @importFrom Rfast Log rep_row Sort
#' @import data.table
#' @importFrom iotools ctapply
#' @importFrom stringi stri_c
#' @export
sequenceStepC1 <- function(amplifiedMolecules, 
                           equalizationAmount, 
                           totalDepth,
                           efficiencyPCR, roundsPCR, efficiencyTag,
                           genes, useUMI) {
  
  
  numCells <- length(amplifiedMolecules)
  amplifiedMoleculesQuant <- equalizeCells(amplifiedMolecules, equalizationAmount)
  amplifiedMolecules <- NULL
  
  taggedMolecules <- lapply(1:numCells, function(x) {
    geneProbs_all <- Rfast::Log(amplifiedMoleculesQuant[[x]] / sum(as.numeric(amplifiedMoleculesQuant[[x]])))
    names(geneProbs_all) <- names(amplifiedMoleculesQuant[[x]])
    totalT <- round(efficiencyTag[x]*sum(as.numeric(amplifiedMoleculesQuant[[x]])))
    tagdM <- try(rmultinom(n=1, size=totalT, prob=exp(geneProbs_all)), silent=T)
    countsT <- as.vector(tagdM)
    names(countsT) <- rownames(tagdM)
    return(countsT)
  })
  amplifiedMoleculesQuant <- NULL
  
  taggedMolecules <- lapply(taggedMolecules, function(x) x * 2) # fragmenting based on C1 guide.
  
  quantifiedMolecules <- amplifyStep(taggedMolecules, genes = genes, efficiencyPCR, roundsPCR, protocol = "C1")
  taggedMolecules <- NULL
  
  
  amplifiedMoleculesQuant_all_list <- lapply(1:numCells, function(x) {
    y = quantifiedMolecules[[x]]
    names(y) <- stringi::stri_c(names(quantifiedMolecules[[x]]), x, sep="__")
    return(y)
  })
  quantifiedMolecules <- NULL
  
  amplifiedMoleculesQuant_all <- do.call(c, amplifiedMoleculesQuant_all_list)
  amplifiedMoleculesQuant_all_list <- NULL
  
  geneProbs_all <- Rfast::Log(amplifiedMoleculesQuant_all / sum(as.numeric(amplifiedMoleculesQuant_all)))
  names(geneProbs_all) <- names(amplifiedMoleculesQuant_all)
  
  maxDepth <- .Machine$integer.max
  if(totalDepth > maxDepth) {
    seqiter <- ceiling(totalDepth / maxDepth)
    counts <- 0
    for (j in 1:(seqiter-1)) {
      counts <- counts + try(rmultinom(n=1, size=maxDepth, prob=exp(geneProbs_all)), silent=T)
    }
    seqDiff <- totalDepth - (maxDepth)*(seqiter-1)
    counts <- counts + try(rmultinom(n=1, size=seqDiff, prob=exp(geneProbs_all)), silent=T)
    
  } else {
    counts <- try(rmultinom(n=1, size=totalDepth, prob=exp(geneProbs_all)), silent=T)
  }
  
  
  count_tab <- NULL
  umi_tab <- NULL
  
  print("Beginning formatting output...")
  
  cnt_split <- data.table::tstrsplit(rownames(counts), split="__", fixed=TRUE)
  cnt_split.df <- data.table(Gene = cnt_split[[1]], 
                             Cell = cnt_split[[2]], Count = counts[,1])
  cnt_split.df$Cell <- as.numeric(cnt_split.df$Cell)
  
  rm(counts)
  rm(cnt_split)
  
  print(paste0("Processing output for ",numCells ," cells."))
  
  cnt_split.df$Cell <- factor(cnt_split.df$Cell)
  srtd_labels <- paste(sort(as.integer(levels(cnt_split.df$Cell))))
  cnt_split.df$Cell <- factor(cnt_split.df$Cell, levels = srtd_labels)
  cnt_split_cell <- split(cnt_split.df, f = cnt_split.df$Cell)
  
  my_tabs <-  lapply(1:length(cnt_split_cell), function(x) {
    X <- cnt_split_cell[[x]]
    X$ugenes <- gsub("@.*","",X$Gene)  
    X <- X[order(X$ugenes),]
    
    count_tab_tempvec <-  iotools::ctapply(X$Count, X$ugenes, sum)
    count_tab_temp <- matrix(count_tab_tempvec)
    rownames(count_tab_temp) <-  names(count_tab_tempvec)
    
    zeroG <- setdiff(genes, rownames(count_tab_temp))
    zeroExpr <- Rfast::rep_row(0, length(zeroG)); rownames(zeroExpr) <- zeroG
    count_tab_temp <- rbind(count_tab_temp, zeroExpr)
    count_tab <- count_tab_temp[Rfast::Sort(rownames(count_tab_temp)),]
    
    if (useUMI==TRUE) {
      x_nonzero <- subset(X, X$Count > 0)
      umi_tab_temptab <- table(x_nonzero$ugenes)
      umi_tab_temp <- matrix(umi_tab_temptab)
      rownames(umi_tab_temp) <- names(umi_tab_temptab)
      zeroG <- setdiff(genes, rownames(umi_tab_temp))
      zeroExpr <- Rfast::rep_row(0, length(zeroG)); rownames(zeroExpr) <- zeroG
      
      umi_tab_temp <- rbind(umi_tab_temp, zeroExpr)
      umi_tab <- umi_tab_temp[Rfast::Sort(rownames(umi_tab_temp)),]
    } else {umi_tab <- NULL}
    
    return(list(count_tab, umi_tab))
  })
  
  count_tab <- do.call(cbind, sapply(my_tabs, function(x) x[1]))
  colnames(count_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
  # rownames(count_tab) <- genes
  
  if (useUMI==TRUE) {
    umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
    colnames(umi_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
    # rownames(umi_tab) <- genes
  }
  
  return(list(counts = count_tab, umi_counts = umi_tab))
}

#' Sequencing step of the Scaffold simulation for 10X/droplet
#' @param capturedMolecules A list of captured molecules for each sample.
#' @param totalDepth The total sequencing depth of the simulated data. If left NULL, this is taken from the \code{sce} object. If more cells are generated than in the original dataset, then the totalDepth will be scaled up accordingly.
#' @param efficiencyPCR A numeric vector representing the efficiency of PCR for each sample.
#' @param roundsPCR An integer indicating the number of PCR amplification rounds.
#' @param efficiencyTag A numeric vector representing the efficiency of tagging for each cell.
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param useUMI A TRUE/FALSE indicating whether the protocol should use UMIs (Unique Molecular Identifiers). Droplet or 10X protocols have this set as TRUE for the default, otherwise FALSE.
#' @importFrom Rfast Log rep_row Sort
#' @import data.table
#' @importFrom iotools ctapply
#' @importFrom stringi stri_c
#' @export
sequenceStep10X <- function(capturedMolecules, totalDepth,
                            efficiencyPCR, roundsPCR, efficiencyTag,
                            genes, useUMI)
{
  print("Rearranging 10X data")
  numCells <- length(capturedMolecules)
  label_capturedMolecules <- lapply(1:length(capturedMolecules), function(x) stringi::stri_c(capturedMolecules[[x]], x, sep="__"))
  
  tab_capturedMolecules <- lapply(label_capturedMolecules, function(y) {
    easyY <- Rfast::rep_row(1, length(y))
    rownames(easyY) <- y
    return(easyY)
  })
  rm(capturedMolecules)
  
  # Now everything is combined and PCR amp:
  longVectorCounts <- unlist(tab_capturedMolecules)
  names(longVectorCounts) <- do.call(c, label_capturedMolecules)
  
  amplifiedMolecules <- amplifyStep(longVectorCounts, genes = genes, efficiencyPCR, roundsPCR, protocol = "10X")
  rm(longVectorCounts)
  
  # Now a fragmentation step happens:
  geneProbs_all <- Rfast::Log(amplifiedMolecules / sum(amplifiedMolecules))
  names(geneProbs_all) <- names(amplifiedMolecules)
  totalT <- round(efficiencyTag*sum(amplifiedMolecules))
  if(totalT > .Machine$integer.max) totalT <- .Machine$integer.max
  tagdM <- try(rmultinom(n=1, size=totalT, prob=exp(geneProbs_all)), silent=T)
  fragments <- as.vector(tagdM)
  names(fragments) <- rownames(tagdM)
  rm(amplifiedMolecules)
  
  # Now we sequence:
  geneProbs_all <- Rfast::Log(fragments / sum(fragments))
  names(geneProbs_all) <- names(fragments)
  maxDepth <- .Machine$integer.max
  if(totalDepth > maxDepth) {
    seqiter <- ceiling(totalDepth / maxDepth)
    counts <- 0
    for (j in 1:(seqiter-1)) {
      counts <- counts + try(rmultinom(n=1, size=maxDepth, prob=exp(geneProbs_all)), silent=T)
    }
    seqDiff <- totalDepth - (maxDepth)*(seqiter-1)
    counts <- counts + try(rmultinom(n=1, size=seqDiff, prob=exp(geneProbs_all)), silent=T)
    
  } else {
    counts <- try(rmultinom(n=1, size=totalDepth, prob=exp(geneProbs_all)), silent=T)
  }
  
  ## Split matrix into nicer output:
  print(paste0("Processing output for ", numCells ," cells."))
  
  cnt_split <- data.table::tstrsplit(rownames(counts), split="__", fixed=TRUE)
  cnt_split.df <- data.table(Gene = cnt_split[[1]], 
                             Cell = cnt_split[[2]], Count = counts[,1])
  cnt_split.df$Cell <- as.numeric(cnt_split.df$Cell)
  
  rm(cnt_split)
  rm(counts)
  
  cnt_split.df$Cell <- factor(cnt_split.df$Cell)
  srtd_labels <- paste(sort(as.integer(levels(cnt_split.df$Cell))))
  cnt_split.df$Cell <- factor(cnt_split.df$Cell, levels = srtd_labels)
  cnt_split_cell <- split(cnt_split.df, f = cnt_split.df$Cell)
  
  
  my_tabs <-  lapply(1:length(cnt_split_cell), function(x) {
    X <- cnt_split_cell[[x]]
    X$ugenes <- gsub("@.*","",X$Gene)  
    X <- X[order(X$ugenes),]
    count_tab_tempvec <-  iotools::ctapply(X$Count, X$ugenes, sum)
    count_tab_temp <- matrix(count_tab_tempvec)
    rownames(count_tab_temp) <-  names(count_tab_tempvec)
    
    zeroG <- setdiff(genes, rownames(count_tab_temp))
    zeroExpr <- Rfast::rep_row(0, length(zeroG)); rownames(zeroExpr) <- zeroG
    count_tab_temp <- rbind(count_tab_temp, zeroExpr)
    count_tab <- count_tab_temp[Rfast::Sort(rownames(count_tab_temp)),]
    
    x_nonzero <- subset(X, X$Count > 0)
    umi_tab_temptab <- table(x_nonzero$ugenes)
    umi_tab_temp <- matrix(umi_tab_temptab)
    rownames(umi_tab_temp) <- names(umi_tab_temptab)
    zeroG <- setdiff(genes, rownames(umi_tab_temp))
    zeroExpr <- Rfast::rep_row(0, length(zeroG)); rownames(zeroExpr) <- zeroG
    
    umi_tab_temp <- rbind(umi_tab_temp, zeroExpr)
    umi_tab <- umi_tab_temp[Rfast::Sort(rownames(umi_tab_temp)),]
    return(list(count_tab, umi_tab))
  })
  
  count_tab <- do.call(cbind, sapply(my_tabs, function(x) x[1]))
  colnames(count_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
  # rownames(count_tab) <- genes
  
  umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
  colnames(umi_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
  # rownames(umi_tab) <- genes
  
  return(list(counts = count_tab, umi_counts = umi_tab))
}