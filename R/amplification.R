# the first amplication
preamplifyStep <- function(capturedMolecules, genes, efficiencyPCR, rounds, typeAMP, useUMI=FALSE){

  lapply(1:length(capturedMolecules), function(x){
	 
	if (useUMI == FALSE) {
		X <- table(capturedMolecules[[x]])}
	else if (useUMI == TRUE) {
		X <- rep(1, length(capturedMolecules[[x]]))
		names(X) <- capturedMolecules[[x]]
	}
	
    if (typeAMP == "PCR") {
      A <- X * (1 + efficiencyPCR[x])^rounds
    }
    if (typeAMP == "IVT") {
      A <- X * (1 + efficiencyPCR[x])*rounds
    }

    zeroG <- setdiff(genes, names(A))
    zeroExpr <- rep(0, length(zeroG)); names(zeroExpr) <- zeroG

    ampMolecules <- c(A, zeroExpr)
    ampMolecules <- ampMolecules[sort(names(ampMolecules))]

    return(ampMolecules)
  })
}
# the second amplification step, this needs to check for protocol and have a separate body for 10x
amplifyStep <- function(capturedMolecules, genes, efficiencyPCR, rounds, protocol){

  if (protocol == "C1")
  {
    lapply(1:length(capturedMolecules), function(x){

      X <- capturedMolecules[[x]]
      A <- X * (1 + efficiencyPCR[x])^rounds

      zeroG <- setdiff(genes, names(A))
      zeroExpr <- rep(0, length(zeroG)); names(zeroExpr) <- zeroG

      ampMolecules <- c(A, zeroExpr)
      ampMolecules <- ampMolecules[sort(names(ampMolecules))]
      return(ampMolecules)
    })
  }
  else if (protocol == "10x" | protocol == "10X") # under construction
  {
    print("Amplification step")

capturedMolecules <- capturedMolecules * (1 + efficiencyPCR)^rounds

    # lapply(1:length(capturedMolecules), function(x){
 # 		easyX <- rep(1, length(capturedMolecules[[x]]))
 # 		names(easyX) <- capturedMolecules[[x]]
 #
 # 	    easyX <- easyX * (1 + efficiencyPCR)^rounds
 #
 # 	    print(paste("Starting cell", x))
 #
 # 		X$ugenes <- gsub("@.*","",X$Gene)
 #        countValues <- X[which(grepl(cellName, names(X)))]
 #        nonZeroNames <- names(countValues)
 #
 #
 #        geneCellNamesForCurrentCell <- paste0(genes, "_", cellName)
 #
 #        zeroGenes <- setdiff(geneCellNamesForCurrentCell, nonZeroNames)
 #        zeroExpr <- rep(0, length(zeroGenes))
 #
 #        names(zeroExpr) <- zeroGenes
 #        countValues <- c(zeroExpr, countValues)
 #        countValues <- countValues[sort(names(countValues))]
 #      return(countValues)
 #    })

  }

}


# equalitzion step
#' @importFrom stats rmultinom
quantCells <- function(amplifiedMolecules, pcntRange) {

  totalM <- round(sapply(amplifiedMolecules, function(x) sum(x)))
  ## Use a negative value here to indicate diluting ALL by same factor...no equalization.
  if (pcntRange < 0) {
    totalM <- totalM / abs(pcntRange)
  } else {
    target <- median(c(min(totalM), quantile(totalM, pcntRange)))
    for(i in 1:length(amplifiedMolecules)) {

      if(totalM[i] > target){
        useMean <- target / totalM[i]
        SF <- rnorm(1, useMean, sd=.01)
      } else {SF = rnorm(1, .95, sd=.01)}
      totalM[i] <- totalM[i] * SF
    }
  }

  inputRange <- totalM
  if (any(inputRange >= 2147483647)) { # largest value in R, just rescaling.
    SCALEALL <- max(inputRange) / 2147483647
    inputRange <- inputRange/SCALEALL
  }

  geneProbs <- lapply(amplifiedMolecules, function(x) log(x) - log(sum(x)))

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
