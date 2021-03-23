#' Estimate scaffold parameters from a SingleCellExperiment.
#'
#' @param sce A SingleCellExperiment object used to estimate the scaffold parameters.
#' @param numCells The number of cells to be used in the simulation. If left NULL, the number of cells in the \code{sce} object is used.
#' @param numGenes The number of genes to be used in the simulation. If left NULL, the number of genes in the \code{sce} object is used.
#' @param geneMeans The mean expression level of each gene. If left NULL, the means are estimated using the \code{sce} object.
#' @param geneTheta The theta parameter used if the negative binomial distribution is chosen to generate the initial counts for each gene. Should be a numeric vector of length \code{numGenes}. If left NULL, this value is estimated from the \code{sce} object.
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param captureEfficiency A value between 0 and 1 indicating the proportion of mRNA molecules successfully captured for the genes in each cell. If left NULL, this value is estimated using the \code{sce} object.
#' @param typeOfAmp The amplification method used in the simulation, defaults to "PCR".
#' @param numFirstAmpCycles The number of cycles in the first amplification stage of the simulation.
#' @param numSecondAmpCycles The number of cycles in the second amplification stage of the simulation.
#' @param firstAmpEfficiency A value between 0 and 1 indicating the efficiency of the first simulated amplification cycle. If set to 1, then all molecules will double each cycle.
#' @param secondAmpEfficiency A value between 0 and 1 indicating the efficiency of the second simulated amplification cycle.
#' @param tagEfficiency A value between 0 and 1 indicating the tagmentation efficiency.
#' @param degree Numeric, determines the amount of perturbation applied to the simulated initial gene counts.
#' @param percentRange A value between 0 and 1 indicating the number of samples that undergo exact dilution in the equalization step of the simulation. A negative value indicates that all samples are diluted by the same factor, resulting in no equalization performed.
#' @param protocol The protocol to model in the simulation.
#' @param totalSD The total sequencing depth of the simulated data. If left NULL, this is taken from the \code{sce} object.
#' @param useUMI A TRUE/FALSE indicating whether the protocol should use UMI (Unique Molecular Identifier).
#' @param model The probability distribution used to generate the initial gene counts. Defaults to 'p' for the Poisson distribution, but can also be set to 'nb' to use the Negative Binomial distribution.
#'
#' @importFrom stats rnorm runif var
#' @importFrom methods new
#' @export
estimateScaffoldParameters <- function(sce,
                                       numCells = NULL,
                                       numGenes = NULL,
                                       geneMeans = NULL,
                                       geneTheta = NULL,
                                       totalTranscripts = NULL,
                                       genes = NULL,
																			 # geneEfficiency = NULL,
                                       captureEfficiency = NULL,
																			 efficiencyRT = NULL,
                                       typeOfAmp = "PCR",
                                       numFirstAmpCycles = 18,
                                       numSecondAmpCycles = 12,
                                       firstAmpEfficiency = NULL,
                                       secondAmpEfficiency = NULL,
                                       tagEfficiency = NULL,
                                       degree = NULL,
                                       percentRange = -1,
                                       protocol = "C1",
                                       totalSD = NULL,
                                       useUMI = FALSE,
                                       model = 'p',
                                       usePops = NULL)
{
  if(is.null(totalTranscripts)) totalTranscripts <- 300000
  # if(is.null(totalTranscripts) & useUMI == FALSE) totalTranscripts <- 300000
  
  # if(is.null(totalTranscripts) & useUMI == TRUE) totalTranscripts <- 10000	
  
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
	# Unlikely to ever use nb, keep for now...
  if (is.null(geneTheta) & model == 'nb') {
    Var <- apply(scaleData, 1, function(x) var(x[x < quantile(x, .9)]))
    Var[is.na(Var)] <- sample(Var[!is.na(Var) & Var > quantile(Var[Var>0], .2, na.rm=T)
                                  & Var < quantile(Var[Var>0], .5, na.rm=T)],
                              sum(is.na(Var)), replace=T)
    Theta <- (geneMeans^2) / (Var - geneMeans)
    Theta[which(Theta <= 0)] <- sample(Theta[Theta > 0 & Theta > quantile(Theta[Theta>0], .1,
                                       na.rm=T) & Theta < quantile(Theta[Theta>0], .5, na.rm=T)],
                                       sum(Theta <= 0, na.rm=T), replace=T)
    geneTheta <- Theta
  } 
  if (is.null(genes)) {
    genes <- rownames(sce)
  }
  if (is.null(captureEfficiency)) {
    captureEfficiency <- -1
  }
	
  if (is.null(firstAmpEfficiency)) {
    firstAmpEfficiency <- Rfast::Rnorm(sum(numCells), .95, .02)
  }
  if (is.null(secondAmpEfficiency)) {
    if (protocol == "C1")
      secondAmpEfficiency <- Rfast::Rnorm(sum(numCells), .95, .02)
    else
      secondAmpEfficiency <- 0.95
  }
  if (is.null(tagEfficiency)) {
    tagEfficiency <- Rfast::Rnorm(sum(numCells), .95, .02)
  }
  
  if(is.null(degree)) {
    degree <- c()
    mycomb <- Rfast::comb_n(sum(numCells), 2)
    
    lib.size <- colSums(counts(sce))
    if (sum(numCells) != length(lib.size)) {
      lib.size <- sample(lib.size, sum(numCells), replace=TRUE)
    }
    
    if (sum(numCells) > 100) {
      mycomb <- mycomb[,sample(1:ncol(mycomb), 100)]
    }
    
    # allComb <- sapply(lib.size, function(x) apply(mycomb,2,function(y) x / lib.size[y[2]]))
    allComb <- apply(mycomb,2,function(y) lib.size[y[1]] / lib.size[y[2]])
    
    degree[1] <- quantile(allComb, .05)	
    degree[2] <- quantile(allComb, .95)
  }
  if (is.null(totalSD)) {
    totalSD <- sum(counts(sce))
  }
  
  if (protocol == "10X") {
    useUMI <- TRUE
  }
  
  return(new("ScaffoldParams",
             numCells = numCells,
             numGenes = numGenes,
             geneMeans = geneMeans,
             geneTheta = geneTheta,
						 totalTranscripts = totalTranscripts,
             genes = genes,
             captureEfficiency = captureEfficiency,
						 efficiencyRT = efficiencyRT,
             typeOfAmp = typeOfAmp,
             numFirstAmpCycles = numFirstAmpCycles,
             numSecondAmpCycles = numSecondAmpCycles,
             firstAmpEfficiency = firstAmpEfficiency,
             secondAmpEfficiency = secondAmpEfficiency,
             tagEfficiency = tagEfficiency,
             degree = degree,
             percentRange = percentRange,
             protocol = protocol,
             totalSD = totalSD,
             useUMI = useUMI,
						 model = model,
             usePops = usePops))
}



# Estimate capture efficiency for simulated data
# method to estimate parameter for capture efficiency
#' @importFrom stats dbinom optimize
#' @export
estimateCaptureEff <- function(Data, compareData) {

  gdetectRate <- rowSums(Data!=0) /ncol(Data)

  Ptest <- rowMeans(Data) / nrow(Data)
  try1 <- split(Rfast::Sort(Ptest), cut(seq_along(Rfast::Sort(Ptest)), 10, labels = FALSE))
  randG <- do.call(c,lapply(1:10, function(x) sample(names(try1[[x]]), 100)))
  minFunc <- function(inGuess){
     tt <- dbinom(0, round(inGuess*nrow(Data)), Ptest[randG], log = TRUE)
     avg.detection.raw = mean(colMeans(compareData!=0))
     X = abs((1-mean(exp(tt), na.rm=T)) - avg.detection.raw)
     return(X)
   }
   simparm <- optimize(minFunc, lower=0, upper=1, tol=1e-10)$minimum
	 getsd <- colMeans(compareData!=0)
	 getsd <- mad(getsd, low=TRUE, constant = 1)
   simparm <- c(simparm, getsd)
   print(simparm) # remove this print later

   esteff <- rnorm(ncol(Data), simparm[1], simparm[2])
   return(esteff)

}




