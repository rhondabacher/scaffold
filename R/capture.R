#' Capture cDNA from cell (cell lysis and conversion of mRNA to cDNA)
#' @param Data A rounded numeric matrix of gene expression counts.
#' @param captureEffCell A vector of values between 0 and 1 to indicate the proportion of mRNA molecules successfully captured for the genes in each cell.
#' @param captureEffGene A vector of values between 0 and 1. If left NULL, all genes are assumed to have equal capture efficiency, and the vector is normalized to sum to 1.
#' @param rtEffCell A vector of values between 0 and 1 to indicate the proportion of mRNA succesfully converted to cDNA.
#' @param rtEffGene A vector with values normalized to sum to 1.
#' @param useUMI A TRUE/FALSE indicating whether the protocol should use UMIs (Unique Molecular Identifiers). Droplet or 10X protocols have this set as TRUE for the default, otherwise FALSE.
#' 
#' @importFrom wrswoR sample_int_rej
#' @export
captureStep <- function(Data, captureEffCell = NULL, captureEffGene = NULL, 
                        rtEffCell = NULL, rtEffGene = NULL, useUMI = FALSE){
  
  if (is.null(captureEffGene)) {
    captureEffGene <- Rfast::rep_col(1, nrow(Data))
    captureEffGene <- as.vector(captureEffGene)
    names(captureEffGene) <- rownames(Data)
  }
  captureEffGene <- captureEffGene / sum(captureEffGene)
  
  if (!is.null(rtEffCell) & is.null(rtEffGene)) {
    rtEffGene <- Rfast::rep_col(1, nrow(Data))
    rtEffGene <- as.vector(rtEffGene)
    names(rtEffGene) <- rownames(Data)
		rtEffGene <- rtEffGene / sum(rtEffGene)
  }
  
  Genes <- rownames(Data)
  
  capturedMolecules <- lapply(seq_len(ncol(Data)), function(x) {
	  
  ## First part is the cell lysis:
	countAll <- Data[,x]
	efficiencyG <- captureEffGene[names(countAll)]
	geneProbs <- rep(efficiencyG, countAll)
	geneProbs <- geneProbs / sum(geneProbs)

	allMolecules <- rep(Genes, countAll)
	totalM <- abs(round(captureEffCell[x]*length(allMolecules)))

	sampledM <- sample_int_rej(n = length(allMolecules), size = totalM, prob = geneProbs)

	sampledM <- allMolecules[sampledM]

  ## Second part is the reverse transcription:
	 if (!is.null(rtEffCell)) {
		countAll <- Rfast::Table(sampledM)
		efficiencyG <- rtEffGene[names(countAll)]
		geneProbs <- rep(efficiencyG, countAll)
		geneProbs <- geneProbs / sum(geneProbs)

		totalM <- abs(round(rtEffCell[x]*length(sampledM)))
		sampledM2 <- sample_int_rej(n = length(sampledM), size = totalM, prob = geneProbs)
		sampledM2 <- sampledM[sampledM2]


	  if (useUMI == TRUE) {
	    outSample <- make.unique(sampledM2, "@")
	  } else {
	    outSample <- sampledM2
	  }

  }
  if (is.null(rtEffCell)){
	
	  if (useUMI == TRUE) {
	    outSample <- make.unique(sampledM, "@")
	  } else {
	    outSample <- sampledM
	  }	
  }
  return(outSample)
 })
 return(capturedMolecules)
}
