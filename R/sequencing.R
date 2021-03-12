sequenceStepC1 <- function(amplifiedMolecules, pcntRange=0, totalSD=50000000,
                         efficiencyPCR, roundsPCR, efficiencyTag=NULL,
						 genes=scaffoldParams@genes, useUMI=FALSE) {


	 amplifiedMoleculesQuant <- quantCells(amplifiedMolecules, pcntRange=pcntRange)
	 amplifiedMolecules <- NULL

	 taggedMolecules <- lapply(1:length(amplifiedMoleculesQuant), function(x) {
	   geneProbs_all <- log(amplifiedMoleculesQuant[[x]] / sum(as.numeric(amplifiedMoleculesQuant[[x]])))
	   totalT <- round(efficiencyTag[x]*sum(as.numeric(amplifiedMoleculesQuant[[x]])))
	   tagdM <- try(rmultinom(n=1, size=totalT, prob=exp(geneProbs_all)), silent=T)
	   countsT <- as.vector(tagdM)
	   names(countsT) <- rownames(tagdM)
	   return(countsT)
	 })
	 amplifiedMoleculesQuant <- NULL

	 taggedMolecules <- lapply(taggedMolecules, function(x) x * 10)

	 quantifiedMolecules <- amplifyStep(taggedMolecules, genes=genes, efficiencyPCR, roundsPCR, protocol = "C1")
	 taggedMolecules <- NULL


	 amplifiedMoleculesQuant_all_list <- lapply(1:length(quantifiedMolecules), function(x) {
	   y = quantifiedMolecules[[x]]
	   names(y) <- paste0(names(quantifiedMolecules[[x]]),"__", x)
	   return(y)
	 })
	 quantifiedMolecules <- NULL

	 amplifiedMoleculesQuant_all <- do.call(c, amplifiedMoleculesQuant_all_list)
	 amplifiedMoleculesQuant_all_list <- NULL

	 geneProbs_all <- log(amplifiedMoleculesQuant_all / sum(as.numeric(amplifiedMoleculesQuant_all)))

	 counts <- try(rmultinom(n=1, size=totalSD, prob=exp(geneProbs_all)), silent=T)

	 count_tab <- NULL
	 umi_tab <- NULL
	 if (useUMI==FALSE) {
	   count_tab <- matrix(counts, nrow=length(genes), ncol=length(efficiencyPCR), byrow=FALSE)
	   rownames(count_tab) <- names(quantifiedMolecules[[1]])
	   colnames(count_tab) <- paste0("Cell_", 1:length(efficiencyPCR))
	 } else if (useUMI==TRUE) {
	   cnt_split <- data.frame(do.call(rbind, strsplit(rownames(counts), split="__")), Count = counts[,1])
	   colnames(cnt_split) <- c("Gene", "Cell", "Count")
	   cnt_split_cell <- split(cnt_split, f = cnt_split$Cell)
	   my_tabs <- lapply(1:length(cnt_split_cell), function(x) {
	     X <- cnt_split_cell[[x]]
	     X$ugenes <- gsub("@.*","",X$Gene)  
	     count_tab <-  tapply(X$Count, X$ugenes, sum)

	     x_nonzero <- subset(X, Count > 0)
	     umi_tab_temp <- table(x_nonzero$ugenes)
	     zeroG <- setdiff(genes, names(umi_tab_temp))
	     zeroExpr <- rep(0, length(zeroG)); names(zeroExpr) <- zeroG

	     umi_tab_temp <- c(umi_tab_temp, zeroExpr)
	     umi_tab <- umi_tab_temp[sort(names(umi_tab_temp))]
	     return(list(count_tab, umi_tab))
	   })

	   count_tab <- do.call(cbind, sapply(my_tabs, function(x) x[1]))
	   colnames(count_tab) <- paste0("Cell_", 1:length(efficiencyPCR))

	   umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
	   colnames(umi_tab) <- paste0("Cell_", 1:length(efficiencyPCR))

	 }

	 return(list(counts = count_tab, umi_counts = umi_tab))


}

# under construction
sequenceStep10X <- function(capturedMolecules, totalSD=50000000,
                            efficiencyPCR, roundsPCR=12, efficiencyTag=NULL,
							genes, useUMI=TRUE)
{
  print("Rearranging 10x data")

  label_capturedMolecules <- lapply(1:length(capturedMolecules), function(x) paste0(capturedMolecules[[x]], "__", x))

  tab_capturedMolecules <- lapply(label_capturedMolecules, function(y) {
    easyY <- rep(1, length(y))
    names(easyY) <- y
    return(easyY)
  })
  
  # Now everything is combined and PCR amp:
  allCounts <- unlist(tab_capturedMolecules)

  amplifiedMolecules <- amplifyStep(longVectorCounts, genes=genes, efficiencyPCR, roundsPCR, protocol = "10X")
  
  # Now a fragmentation step happens:
  geneProbs_all <- log(amplifiedMolecules / sum(amplifiedMolecules))
  totalT <- round(efficiencyTag*sum(amplifiedMolecules))
  tagdM <- try(rmultinom(n=1, size=totalT, prob=exp(geneProbs_all)), silent=T)
  fragments <- as.vector(tagdM)
  names(fragments) <- rownames(tagdM)

  # Now we sequence:
  geneProbs_all <- log(fragments / sum(fragments))
  counts <- try(rmultinom(n=1, size=totalSD, prob=exp(geneProbs_all)), silent=T)
  
  ## Split matrix into nicer output:

  cnt_split <- data.frame(do.call(rbind, strsplit(rownames(counts), split="__")), Count = counts[,1])
  colnames(cnt_split) <- c("Gene", "Cell", "Count")
  cnt_split$Cell <- factor(cnt_split$Cell)
  cnt_split$Cell <- factor(cnt_split$Cell, levels = order(levels(cnt_split$Cell)))
  cnt_split_cell <- split(cnt_split, f = cnt_split$Cell)
  my_tabs <- lapply(1:length(cnt_split_cell), function(x) {
    X <- cnt_split_cell[[x]]
    X$ugenes <- gsub("@.*","",X$Gene)  
    count_tab_temp <-  tapply(X$Count, X$ugenes, sum)
    zeroG <- setdiff(genes, names(count_tab_temp))
    zeroExpr <- rep(0, length(zeroG)); names(zeroExpr) <- zeroG
    count_tab_temp <- c(count_tab_temp, zeroExpr)
    count_tab <- count_tab_temp[sort(names(count_tab_temp))]

    x_nonzero <- subset(X, Count > 0)
    umi_tab_temp <- table(x_nonzero$ugenes)
    zeroG <- setdiff(genes, names(umi_tab_temp))
    zeroExpr <- rep(0, length(zeroG)); names(zeroExpr) <- zeroG

    umi_tab_temp <- c(umi_tab_temp, zeroExpr)
    umi_tab <- umi_tab_temp[sort(names(umi_tab_temp))]
    return(list(count_tab, umi_tab))
  })

  count_tab <- do.call(cbind, sapply(my_tabs, function(x) x[1]))
  colnames(count_tab) <- paste0("Cell_", 1:length(capturedMolecules))

  umi_tab <- do.call(cbind, sapply(my_tabs, function(x) x[2]))
  colnames(umi_tab) <- paste0("Cell_", 1:length(capturedMolecules))

 return(list(counts = count_tab, umi_counts = umi_tab))
}

