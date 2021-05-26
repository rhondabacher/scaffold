#' Sequencing step of the Scaffold simulation
#' @importFrom Rfast Log rep_row Sort
#' @import data.table
#' @importFrom iotools ctapply
#' @importFrom stringi stri_c
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
	 
	 counts <- try(rmultinom(n=1, size=totalDepth, prob=exp(geneProbs_all)), silent=T)
	
	 count_tab <- NULL
	 umi_tab <- NULL
	 
	 print("Beginning formatting output...")
	 
	 cnt_split <- unlist(strsplit(rownames(counts), split="__", fixed=TRUE))
	 cnt_split <- data.table(Gene = (cnt_split)[c(TRUE,FALSE)], Cell = (cnt_split)[c(FALSE,TRUE)], Count = counts[,1])

	 print(paste0("Processing output for ",numCells ," cells."))
	 cnt_split$Cell <- factor(cnt_split$Cell)
	 cnt_split$Cell <- factor(cnt_split$Cell, levels = order(levels(cnt_split$Cell)))
	 cnt_split_cell <- split(cnt_split, f = cnt_split$Cell)
   
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
	   rownames(count_tab) <- genes
	   
		 if (useUMI==TRUE) {
			 umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
			 colnames(umi_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
			 rownames(umi_tab) <- genes
		 }
	  
	 return(list(counts = count_tab, umi_counts = umi_tab))
}

#' Sequencing step of the Scaffold simulation for 10X/droplet
#' @importFrom Rfast Log rep_row Sort
#' @import data.table
#' @importFrom iotools ctapply
#' @importFrom stringi stri_c
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
  counts <- try(rmultinom(n=1, size=totalDepth, prob=exp(geneProbs_all)), silent=T)
  
  ## Split matrix into nicer output:
  print(paste0("Processing output for ", numCells ," cells."))
	cnt_split <- unlist(strsplit(rownames(counts), split="__", fixed=TRUE))
	cnt_split <- data.table(Gene = (cnt_split)[c(TRUE,FALSE)], Cell = (cnt_split)[c(FALSE,TRUE)], Count = counts[,1])
  cnt_split$Cell <- factor(cnt_split$Cell)
  cnt_split$Cell <- factor(cnt_split$Cell, levels = order(levels(cnt_split$Cell)))
  cnt_split_cell <- split(cnt_split, f = cnt_split$Cell)
	

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
  rownames(count_tab) <- genes
  
  umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
  colnames(umi_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
  rownames(umi_tab) <- genes
  return(list(counts = count_tab, umi_counts = umi_tab))
}