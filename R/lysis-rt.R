                                        # this is step 2
#' @importFrom wrswoR sample_int_rej
captureStep <- function(Data, captureEffCell = NULL, captureEffGene = NULL, 
							rtEffCell = NULL, rtEffGene = NULL, useUMI = FALSE){
  
  if (is.null(captureEffGene) == TRUE) {
    captureEffGene <- rep(1, nrow(Data))
  }
  captureEffGene <- captureEffGene / sum(captureEffGene)
  names(captureEffGene) <- rownames(Data)
  if (is.null(rtEffGene) == TRUE) {
    rtEffGene <- rep(1, nrow(Data))
  }
  rtEffGene <- rtEffGene / sum(rtEffGene)
  names(rtEffGene) <- rownames(Data)
  
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
		countAll <- table(sampledM)
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






#                                         # this is step 2
# #' @importFrom wrswoR sample_int_rej
# lysisStep <- function(Data, efficiencyCell = NULL, efficiencyGene = NULL, useUMI=FALSE){
#
#   if (is.null(efficiencyGene) == TRUE) {
#     efficiencyGene <- rep(1, nrow(Data))
#   }
#
#   efficiencyGene <- efficiencyGene / sum(efficiencyGene)
#   names(efficiencyGene) <- rownames(Data)
#
#   Genes <- rownames(Data)
#
#   capturedMolecules <- lapply(seq_len(ncol(Data)), function(x) {
#     countAll <- Data[,x]
#     efficiencyG <- efficiencyGene[names(countAll)]
#     geneProbs <- rep(efficiencyG, countAll)
#     geneProbs <- geneProbs / sum(geneProbs)
#
#     allMolecules <- rep(Genes, countAll)
#     totalM <- abs(round(efficiencyCell[x]*length(allMolecules)))
#
#     sampledM <- sample_int_rej(n = length(allMolecules), size = totalM, prob = geneProbs)
#
#     if (useUMI == TRUE) {
#       sampledM <- make.unique(allMolecules[sampledM])
#     } else {
#       sampledM <- allMolecules[sampledM]
#     }
#
#     return(sampledM)
#   })
#   return(capturedMolecules)
# }
#
# revtStep <- function(Data, efficiencyCell = NULL, efficiencyGene = NULL,
#                      Genes=NULL){
#
#   if(is.null(efficiencyGene) == TRUE) {
#     efficiencyGene <- rep(1, length(Genes))
#   }
#   efficiencyGene <- efficiencyGene / sum(efficiencyGene)
#   names(efficiencyGene) <- Genes
#
#   capturedMolecules <- lapply(1:length(Data), function(x) {
#     countAll <- table(Data[[x]])
#
#     efficiencyG <- efficiencyGene[names(countAll)]
#     geneProbs <- rep(efficiencyG, countAll)
#     geneProbs <- geneProbs / sum(geneProbs)
#
#     totalM <- abs(round(efficiencyCell[x]*length(Data[[x]])))
#     sampledM <- sample_int_rej(n = length(Data[[x]]), size = totalM, prob = geneProbs)
#     sampledM <- Data[[x]][sampledM]
#
#     return(sampledM)
#   })
# }
