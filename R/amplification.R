# the first amplication
#' @export
preamplifyStep <- function(capturedMolecules, genes, efficiencyPCR, rounds, typeAMP, useUMI=FALSE){

  lapply(1:length(capturedMolecules), function(x){
	 
	if (useUMI == FALSE) {
		X <- Rfast::Table(capturedMolecules[[x]])}
	else if (useUMI == TRUE) {
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

# the second amplification step, this needs to check for protocol and have a separate body for 10x
#' @export
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
  else if (protocol == "10x" | protocol == "10X") # under construction
  {
    print("Amplification step")

	capturedMolecules <- capturedMolecules * (1 + efficiencyPCR)^rounds

  }

}


# equalitzion step
#' @importFrom stats rmultinom
#' @export
quantCells <- function(amplifiedMolecules, pcntRange) {

  totalM <- round(sapply(amplifiedMolecules, function(x) sum(x)))
  ## Use a negative value here to indicate diluting ALL by same factor...no equalization.
  if (pcntRange < 0) {
    totalM <- totalM * rnorm(length(totalM), .95, sd=.01)
  } else {
    target <- median(c(min(totalM), quantile(totalM, pcntRange)))
    for(i in 1:length(amplifiedMolecules)) {

      if(totalM[i] > target){
        useMean <- target / totalM[i]
        SF <- abs(rnorm(1, useMean, sd=.01))
      } else {SF = rnorm(1, .95, sd=.01)}
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
