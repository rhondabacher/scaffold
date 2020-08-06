#' Estimate scaffold parameters from a SingleCellExperiment.
#'
#' @param sce A SingleCellExperiment object used to estimate the scaffold parameters.
#' @param numCells The number of cells to be used in the simulation. If left NULL, the number of cells in the \code{sce} object is used.
#' @param numGenes The number of genes to be used in the simulation. If left NULL, the number of genes in the \code{sce} object is used.
#' @param geneMeans The mean expression level of each gene. If left NULL, the means are estimated using the \code{sce} object.
#' @param geneTheta XXTODOXX
#' @param genes A vector of names for each gene. If left NULL, the gene names from the \code{sce} object are used.
#' @param captureEfficiency XXTODOXX
#' @param typeOfAmp The amplification method used in the simulation, defaults to "PCR".
#' @param numFirstAmpCycles The number of cycles in the first amplification stage of the simulation.
#' @param numSecondAmpCycles The number of cycles in the second amplification stage of the simulation.
#' @param firstAmpEfficiency A value between 0 and 1 indicating the efficiency of the first simulated amplification cycle.
#' @param secondAmpEfficiency A value between 0 and 1 indicating the efficiency of the second simulated amplification cycle.
#' @param tagEfficiency XXTODOXX
#' @param degree Numeric, determines the amount of perturbation applied to the simulated initial gene counts.
#' @param percentRange A value between 0 and 1 indicating the number of samples that undergo exct diluation in the equalization step of the simulation. A negative value indicates that all samples are diluated by the same factor, resulting in no equalization performed.
#' @param protocol The protocol to model in the simulation.
#' @param totalSD The total sequencing depth of the simulated data. If left NULL, this is taken from the \code{sce} object.
#'
#' @importFrom stats rnorm runif var
#' @importFrom methods new
#' @export
estimateScaffoldParameters <- function(sce,
                                     numCells = NULL,
                                     numGenes = NULL,
                                     geneMeans = NULL,
                                     geneTheta = NULL,
                                     genes = NULL,
                                     captureEfficiency = NULL,
                                     typeOfAmp = "PCR",
                                     numFirstAmpCycles = 18,
                                     numSecondAmpCycles = 12,
                                     firstAmpEfficiency = NULL,
                                     secondAmpEfficiency = NULL,
                                     tagEfficiency = NULL,
                                     degree = 5,
                                     percentRange = -1,
                                     protocol = "C1",
                                     totalSD = NULL)
{
  SF <- colSums(counts(sce)) / 500000
  NORMTRY <- t(t(counts(sce)) / SF)
  if (is.null(numCells))
    numCells <- ncol(counts(sce))
  if (is.null(numGenes))
    numGenes <- nrow(counts(sce))
  if (is.null(geneMeans))
  {
    Mu <- apply(NORMTRY, 1, function(x) mean(x))
    Mu[is.na(Mu)] <- min(Mu, na.rm=T)
    geneMeans <- Mu
  }
  if (is.null(geneTheta))
  {
    Var <- apply(NORMTRY, 1, function(x) var(x[x < quantile(x, .9)]))
    Var[is.na(Var)] <- sample(Var[!is.na(Var) & Var > quantile(Var[Var>0], .2, na.rm=T)
                                  & Var < quantile(Var[Var>0], .5, na.rm=T)],
                              sum(is.na(Var)), replace=T)
    Theta <- (geneMeans^2) / (Var - geneMeans)
    Theta[which(Theta <= 0)] <- sample(Theta[Theta > 0 & Theta > quantile(Theta[Theta>0], .1,
                                                                          na.rm=T) & Theta < quantile(Theta[Theta>0], .5, na.rm=T)],
                                       sum(Theta <= 0, na.rm=T), replace=T)
    geneTheta <- Theta
  }
  if (is.null(genes))
  {
    genes <- rownames(sce)
  }
  if (is.null(captureEfficiency))
  {
    captureEfficiency <- -1
  }
  if (is.null(firstAmpEfficiency))
  {
    firstAmpEfficiency <- rnorm(numCells, .90, .02)
  }
  if (is.null(secondAmpEfficiency))
  {
    if (protocol == "C1")
      secondAmpEfficiency <- rnorm(numCells, .90, .02)
    else
      secondAmpEfficiency <- 0.9
  }
  if (is.null(tagEfficiency))
  {
    tagEfficiency <- runif(numCells, .95, 1)
  }
  if (is.null(totalSD))
  {
    totalSD <- sum(counts(sce))
  }
  return(new("ScaffoldParams",
             numCells = numCells,
             numGenes = numGenes,
             geneMeans = geneMeans,
             geneTheta = geneTheta,
             genes = genes,
             captureEfficiency = captureEfficiency,
             typeOfAmp = typeOfAmp,
             numFirstAmpCycles = numFirstAmpCycles,
             numSecondAmpCycles = numSecondAmpCycles,
             firstAmpEfficiency = firstAmpEfficiency,
             secondAmpEfficiency = secondAmpEfficiency,
             tagEfficiency = tagEfficiency,
             degree = degree,
             percentRange = percentRange,
             protocol = protocol,
             totalSD = totalSD))
}
