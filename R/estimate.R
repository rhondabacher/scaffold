#' Estimate scaffold parameters from a SingleCellExperiment.
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
