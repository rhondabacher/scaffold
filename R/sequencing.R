#' @export
sequenceStepC1 <- function(amplifiedMolecules, pcntRange=0, totalSD=50000000,
                          efficiencyPCR, roundsPCR, efficiencyTag=NULL,
						              genes=scaffoldParams@genes, 
                                      useUMI=FALSE,
                                      cores=1) {


   numCells <- length(amplifiedMolecules)
	 amplifiedMoleculesQuant <- quantCells(amplifiedMolecules, pcntRange=pcntRange)
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

	 taggedMolecules <- lapply(taggedMolecules, function(x) x * 2) # based on C1 guide.

	 quantifiedMolecules <- amplifyStep(taggedMolecules, genes=genes, efficiencyPCR, roundsPCR, protocol = "C1")
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
	 
	 counts <- try(rmultinom(n=1, size=totalSD, prob=exp(geneProbs_all)), silent=T)
	
	 count_tab <- NULL
	 umi_tab <- NULL
	 
	 print("Beginning formatting output...")
	 
	 library(data.table)
	 cnt_split <- unlist(strsplit(rownames(counts), split="__", fixed=TRUE))
	 cnt_split <- data.table(Gene = (cnt_split)[c(TRUE,FALSE)], Cell = (cnt_split)[c(FALSE,TRUE)], Count = counts[,1])

	 print(paste0("Processing output for ",numCells ," cells."))
	 cnt_split$Cell <- factor(cnt_split$Cell)
	 cnt_split$Cell <- factor(cnt_split$Cell, levels = order(levels(cnt_split$Cell)))
	 cnt_split_cell <- split(cnt_split, f = cnt_split$Cell)
   
     cl <- parallel::makeForkCluster(cores)
     doParallel::registerDoParallel(cl)
     ind.block <- bigstatsr:::CutBySize(length(cnt_split_cell), nb = cores)

     library(foreach)

     my_tabs <- foreach(i = rows_along(ind.block)) %dopar% {
       vals <- bigstatsr:::seq2(ind.block[i, ])
       lapply(vals, function(x) {
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
				     x_nonzero <- subset(X, Count > 0)
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
	} 
       my_tabs <- if(length(my_tabs) == cores) { do.call(c, my_tabs) }
     
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

# under construction
#' @export
sequenceStep10X <- function(capturedMolecules, totalSD=50000000,
                            efficiencyPCR, roundsPCR=12, efficiencyTag,
                            genes, useUMI=TRUE, cores = 1)
{
  print("Rearranging 10x data")
  numCells <- length(capturedMolecules)
  label_capturedMolecules <- lapply(1:length(capturedMolecules), function(x) stringi::stri_c(capturedMolecules[[x]], x, sep="__"))
  
  tab_capturedMolecules <- lapply(label_capturedMolecules, function(y) {
    easyY <- Rfast::rep_row(1, length(y))
    rownames(easyY) <- y
    return(easyY)
  })
  
  # Now everything is combined and PCR amp:
  longVectorCounts <- unlist(tab_capturedMolecules)
  names(longVectorCounts) <- do.call(c, label_capturedMolecules)
  
  amplifiedMolecules <- amplifyStep(longVectorCounts, genes=genes, efficiencyPCR, roundsPCR, protocol = "10X")
  
	print("Beginning formatting output...")
  # Now a fragmentation step happens:
  geneProbs_all <- Rfast::Log(amplifiedMolecules / sum(amplifiedMolecules))
  names(geneProbs_all) <- names(amplifiedMolecules)
  totalT <- round(efficiencyTag*sum(amplifiedMolecules))
  if(totalT > .Machine$integer.max) totalT <- .Machine$integer.max
  tagdM <- try(rmultinom(n=1, size=totalT, prob=exp(geneProbs_all)), silent=T)
  fragments <- as.vector(tagdM)
  names(fragments) <- rownames(tagdM)
  
  # Now we sequence:
  geneProbs_all <- Rfast::Log(fragments / sum(fragments))
  names(geneProbs_all) <- names(fragments)
  counts <- try(rmultinom(n=1, size=totalSD, prob=exp(geneProbs_all)), silent=T)
  
  ## Split matrix into nicer output:
  print(paste0("Processing output for ",numCells ," cells."))
	cnt_split <- unlist(strsplit(rownames(counts), split="__", fixed=TRUE))
	cnt_split <- data.table(Gene = (cnt_split)[c(TRUE,FALSE)], Cell = (cnt_split)[c(FALSE,TRUE)], Count = counts[,1])
  cnt_split$Cell <- factor(cnt_split$Cell)
  cnt_split$Cell <- factor(cnt_split$Cell, levels = order(levels(cnt_split$Cell)))
  cnt_split_cell <- split(cnt_split, f = cnt_split$Cell)
	
  cl <- parallel::makeForkCluster(cores)
  doParallel::registerDoParallel(cl)
  ind.block <- bigstatsr:::CutBySize(length(cnt_split_cell), nb = cores)

  library(foreach)

  my_tabs <- foreach(i = rows_along(ind.block)) %dopar% {
    vals <- bigstatsr:::seq2(ind.block[i, ])
    lapply(vals, function(x) {
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
  
	  x_nonzero <- subset(X, Count > 0)
	  umi_tab_temptab <- table(x_nonzero$ugenes)
	  umi_tab_temp <- matrix(umi_tab_temptab)
	  rownames(umi_tab_temp) <- names(umi_tab_temptab)
	  zeroG <- setdiff(genes, rownames(umi_tab_temp))
	  zeroExpr <- Rfast::rep_row(0, length(zeroG)); rownames(zeroExpr) <- zeroG
  
	  umi_tab_temp <- rbind(umi_tab_temp, zeroExpr)
	  umi_tab <- umi_tab_temp[Rfast::Sort(rownames(umi_tab_temp)),]
	  return(list(count_tab, umi_tab))
	})
  }
  my_tabs <- if(length(my_tabs) == cores) { do.call(c, my_tabs) }
  count_tab <- do.call(cbind, sapply(my_tabs, function(x) x[1]))
  colnames(count_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
  rownames(count_tab) <- genes
  
  umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
  colnames(umi_tab) <- stringi::stri_c("Cell", 1:numCells, sep="_")
  rownames(umi_tab) <- genes
  return(list(counts = count_tab, umi_counts = umi_tab))
}