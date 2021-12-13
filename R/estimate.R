#' Estimate or set all parameters for Scaffold simulation.
#'
#' @param sce A SingleCellExperiment object used to estimate some parameters for Scaffold. A data matrix is also acceptable as input. This can be left as NULL if all other parameters are passed in.
#' @param sceUMI Whether or not the data in the SCE object are UMI counts.
#' @param numCells The number of cells to be used in the simulation. If left NULL, the number of cells in the \code{sce} object is used. If simulating multiple populations, this should be a vector.
#' @param numGenes The number of genes to be used in the simulation. If left NULL, the number of genes in the \code{sce} object is used.
#' @param geneMeans The mean expression level of each gene. If left NULL, the means are estimated from the \code{sce} object.
#' @param totalTranscripts The total number of transcripts per cell. If left NULL, the deafult is 300,000.
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param protocol The protocol to model in the simulation (accepted input is: C1, droplet, 10X).
#' @param useUMI A TRUE/FALSE indicating whether the protocol should use UMIs (Unique Molecular Identifiers). Droplet or 10X protocols have this set as TRUE for the default, otherwise FALSE.
#' @param popHet a vector of length two to indicate the lower and upper bounds of the amount of perturbation applied to the simulated initial gene counts. This represents natural population heterogeneity. If NULL, scaffold will estimate from the data. To simulate a homogenous population, set popHet=c(1,1).
#' @param geneEfficiency This parameter is not currently used, but may be implemented in a future release. 
#' @param captureEfficiency A vector of values between 0 and 1 to indicate the proportion of mRNA molecules successfully captured for the genes in each cell. If left NULL, this value is estimated using the \code{sce} object.
#' @param efficiencyRT If left NULL (default), this step of the protocol is skipped. Otherwise, the user can specify a vector of values between 0 and 1 to indicate the proportion of mRNA succesfully converted to cDNA.
#' @param typeOfAmp The amplification method used in the simulation, defaults to "PCR", "IVT" is another accepted value.
#' @param numPreAmpCycles The number of cycles to use in the pre-amplification or first amplification stage of the simulation.
#' @param numAmpCycles The number of cycles to use in the second amplification stage of the simulation, or only amplification for droplet and 10X protocols.
#' @param preAmpEfficiency A vector of values between 0 and 1 indicating the efficiency of the first simulated amplification cycle. If set to 1, then all molecules will double each cycle.
#' @param ampEfficiency A vector of values between 0 and 1 indicating the efficiency of the second simulated amplification cycle. For droplet/10X protocols this should be a vector of length one (a single value).
#' @param tagEfficiency A value between 0 and 1 indicating the tagmentation efficiency.
#' @param equalizationAmount A value between 0 and 1 indicating the q* to determine the number of samples that undergo dilution in the equalization step of the simulation. A value of 0 indicates all cells are diluted to the smallest concentration and a value of 1 indicates no equalization is performed.
#' @param totalDepth The total sequencing depth of the simulated data. If left NULL, this is taken from the \code{sce} object.
#' @param usePops This should be a named list with elements: propGenes, fc_mean, fc_sd. The elements are vectors with length one less than the number of cell populations. propGenes indicates the proportion of genes having distinct expression compared to the first cell population. fc_mean and fc_sd control each populations fold-change mean and standard deviation.
#' @param useDynamic This should be a named list with elements: propGenes, fc_mean, fc_sd. The elements are vectors with length one less than the number of cell populations. propGenes indicates the proportion of genes having distinct expression compared to the first cell population. fc_mean and fc_sd control each populations fold-change mean and standard deviation.
#'
#' @importFrom stats rnorm runif var
#' @importFrom methods new
#' @importFrom Rfast Rnorm comb_n
#' @export
estimateScaffoldParameters <- function(sce = NULL, sceUMI = FALSE, numCells = NULL, numGenes = NULL, geneMeans = NULL,
      totalTranscripts = NULL, genes = NULL, protocol = "C1", useUMI = FALSE, 
			popHet = NULL,
			geneEfficiency = NULL, 
			captureEfficiency = NULL, efficiencyRT = NULL, typeOfAmp = "PCR", numPreAmpCycles = 18,
			numAmpCycles = 12, preAmpEfficiency = NULL, ampEfficiency = NULL, tagEfficiency = NULL,
			equalizationAmount = 1, totalDepth = NULL,
      usePops = NULL, useDynamic = NULL)
{
	
	## Setting up the default values:
	
  if (is.null(totalTranscripts)) totalTranscripts <- 300000

		# Divide out depth before estimating the means.
  scaleFactor <- colSums(counts(sce)) / totalTranscripts
  scaleData <- t(t(counts(sce)) / scaleFactor)
  
  if (is.null(numCells)) {
    numCells <- ncol(counts(sce))
  }  
  if (is.null(numGenes)) {
    numGenes <- nrow(counts(sce))
  }
  if (is.null(geneMeans)) {
    Mu <- rowMeans(scaleData)
    Mu[is.na(Mu)] <- min(Mu, na.rm=T)
    geneMeans <- Mu
  }
  if (is.null(genes)) {
    genes <- rownames(sce)
  }
	if (protocol=="droplet") protocol <- "10X"

  if (is.null(equalizationAmount)) {
    equalizationAmount <- 1 # No equalization is default
  }
  if (is.null(preAmpEfficiency)) {
    preAmpEfficiency <- Rfast::Rnorm(sum(numCells), .95, .02)
  }
  if (is.null(ampEfficiency)) {
    if (protocol == "C1")
      ampEfficiency <- Rfast::Rnorm(sum(numCells), .95, .02)
    else
      ampEfficiency <- 0.95
  }
  if (is.null(tagEfficiency)) {
		if (protocol == "C1")
    tagEfficiency <- Rfast::Rnorm(sum(numCells), .95, .02)
		else
			tagEfficiency <- .95
  }
  
  if (is.null(popHet)) {
    popHet <- c()
    mycomb <- Rfast::comb_n(sum(numCells), 2)
    lib_size <- colSums(counts(sce))
    if (sum(numCells) != length(lib_size)) {
      lib_size <- sample(lib_size, sum(numCells), replace=TRUE)
    }
    
    if (sum(numCells) > 100) {
      mycomb <- mycomb[,sample(seq(1,ncol(mycomb)), 100)]
    }
  
    allComb <- apply(mycomb,2,function(y) lib_size[y[1]] / lib_size[y[2]])
    
    popHet[1] <- quantile(allComb, .05)	
    popHet[2] <- quantile(allComb, .95)
  }
  if (is.null(totalDepth)) {
    totalDepth <- sum(counts(sce))
  }
	
  if (protocol %in% c("droplet","10X")) {
    useUMI <- TRUE
  }
  
  return(new("ScaffoldParams",
						 sceUMI = sceUMI,
             numCells = numCells,
             numGenes = numGenes,
             geneMeans = geneMeans,
						 totalTranscripts = totalTranscripts,
						 genes = genes,
             protocol = protocol,
						 useUMI = useUMI,
						 popHet = popHet,
						 geneEfficiency = geneEfficiency,
						 captureEfficiency = captureEfficiency,
						 efficiencyRT = efficiencyRT,
             typeOfAmp = typeOfAmp,		 
             numPreAmpCycles = numPreAmpCycles,
             numAmpCycles = numAmpCycles,
             preAmpEfficiency = preAmpEfficiency,
             ampEfficiency = ampEfficiency,
             tagEfficiency = tagEfficiency,
             equalizationAmount = equalizationAmount,
             totalDepth = totalDepth,
             usePops = usePops,
						 useDynamic = useDynamic))
}



#' Estimate capture efficiency for simulated data
#' method to estimate parameter for capture efficiency
#' @param Data This is the original data from the sce.
#' @param compareData This is the initial mRNA counts.
#' @param protocol Which scRNA-seq protocol is being simulated.
#' @param fromUMI whether or not the original data are UMI counts.
#' @import stats
estimateCaptureEff <- function(Data, compareData, protocol, fromUMI) {

  gdetectRate <- rowSums(Data!=0) /ncol(Data)

  Pweight <- rowMeans(Data) / nrow(Data)
  splitG <- split(Rfast::Sort(Pweight), cut(seq_along(Rfast::Sort(Pweight)), 10, labels = FALSE))
  randG <- do.call(c,lapply(1:10, function(x) sample(names(splitG[[x]]), 100)))

   if (protocol=="10X" | fromUMI == TRUE) {
       minFuncUMI <- function(inGuess){
         tt <- pbinom(1, round(inGuess*nrow(Data)), Pweight[randG], log.p = TRUE)
         avg.detection.raw = mean(colMeans(compareData > 0))
         X = abs((1 - mean(exp(tt), na.rm=T)) - avg.detection.raw)
         return(X)
       }
	  simparm <- optimize(minFuncUMI, lower=0, upper=1, tol=1e-10)$minimum
	 } else if (protocol=="C1") {
         minFuncC1 <- function(inGuess){
            tt <- dbinom(0, round(inGuess*nrow(Data)), Pweight[randG], log = TRUE)
            avg.detection.raw = mean(colMeans(compareData > 0))
            X = abs((1 - mean(exp(tt), na.rm=T)) - avg.detection.raw)
            return(X)
          }
         simparm <- optimize(minFuncC1, lower=0, upper=1, tol=1e-10)$minimum
     }
	 getsd <- colMeans(compareData != 0)
	 getsd <- mad(getsd, low = TRUE, constant = 1)
   simparm <- c(simparm, getsd)

   esteff <- abs(rnorm(ncol(Data), simparm[1], simparm[2]))
   return(esteff)

}




