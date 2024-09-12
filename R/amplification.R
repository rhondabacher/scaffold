#' Pre-amplification step
#' @param capturedMolecules A truncated lists of genes.
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param efficiencyPCR A numeric vector representing the efficiency of PCR for each sample.
#' @param rounds An integer specifying the number of PCR or IVT amplification rounds to perform.
#' @param typeAMP The amplification method used in the simulation, defaults to "PCR", "IVT" is another accepted value.
#' @param useUMI A TRUE/FALSE indicating whether the protocol should use UMIs (Unique Molecular Identifiers). Droplet or 10X protocols have this set as TRUE for the default, otherwise FALSE.
#' 
#' @importFrom Rfast Table rep_col
preamplifyStep <- function(capturedMolecules, genes, efficiencyPCR, rounds, typeAMP, useUMI){

  lapply(1:length(capturedMolecules), function(x){
	 
		if (useUMI == FALSE) {
			X <- Rfast::Table(capturedMolecules[[x]])
		} else if (useUMI == TRUE) {
			X <- Rfast::rep_col(1, length(capturedMolecules[[x]]))
			X <- as.vector(X)
			names(X) <- capturedMolecules[[x]]
		}
	  if (typeAMP == "PCR") {
	      A <- X * (1 + efficiencyPCR[x])^rounds
	  }
	  if (typeAMP == "IVT") {
	     A <- X * (1 + efficiencyPCR[x])*rounds
	  }

		ampMolecules <- A
	return(ampMolecules)
	})
}

#' The second amplification step
#' @param capturedMolecules A truncated lists of genes.
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param efficiencyPCR A numeric vector representing the efficiency of PCR for each sample.
#' @param rounds An integer specifying the number of PCR or IVT amplification rounds to perform.
#' @param protocol The protocol to model in the simulation (accepted input is: C1, droplet, 10X).
amplifyStep <- function(capturedMolecules, genes, efficiencyPCR, rounds, protocol){

  if (protocol == "C1")
  {
    lapply(1:length(capturedMolecules), function(x){

      X <- capturedMolecules[[x]]
      A <- X * (1 + efficiencyPCR[x])^rounds

      ampMolecules <- A 
      return(ampMolecules)
    })
  }
  else if (protocol == "10x" | protocol == "10X") 
  {

	capturedMolecules <- capturedMolecules * (1 + efficiencyPCR)^rounds

  }

}


# equalitzion step
#' @importFrom stats rmultinom
equalizeCells <- function(amplifiedMolecules, equalizationAmount) {

  totalM <- round(sapply(amplifiedMolecules, function(x) sum(x)))

  if (equalizationAmount == 1) { # no equalization
    totalM <- totalM * rnorm(length(totalM), .95, sd=.01)
  } else {
    target <- median(c(min(totalM), quantile(totalM, equalizationAmount)))
    for(i in 1:length(amplifiedMolecules)) {

      if(totalM[i] > target){
        useMean <- target / totalM[i]
        SF <- abs(rnorm(1, useMean, sd = .01))
      } else {SF = rnorm(1, .95, sd = .01)}
			SF[SF > 1] <- 1
      totalM[i] <- totalM[i] * SF
    }
  }
  inputRange <- totalM
  if (any(inputRange >= .Machine$integer.max)) { # largest value in R, just rescaling.
    SCALEALL <- max(inputRange) / .Machine$integer.max
    inputRange <- inputRange/SCALEALL
  }
  geneProbs <- lapply(amplifiedMolecules, function(x) log(x/sum(x)))

  countsList <- lapply(1:length(geneProbs), function(x) {
    geneProbsUse = geneProbs[[x]]
    counts <- try(rmultinom(n=1, size=inputRange[x], prob=exp(geneProbsUse)), silent=T)
    return(counts)
  })
  
  countsList <- lapply(countsList, function(x) if(class(x)[1] != "try-error") {
    countsX <- as.vector(x)
    names(countsX) <- rownames(x)
    return(countsX)})

  return(countsList)
}
