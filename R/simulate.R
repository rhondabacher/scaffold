#' Run a scaffold simulation.
#'
#' @param scaffoldParams An object of class ScaffoldParams.
#' @param originalSCE The SingleCellExperiment used to create the scaffoldParams object.
#' @param model The probability distribution used to generate the initial gene counts. Defaults to 'p' for the Poisson distribution, but can also be set to 'nb' to use the Negative Binomial distribution.
#' @param inputInitial A optional matrix of initial gene counts which should have the same dimension indicated in the \code{scaffoldParams} parameter. If left NULL, the initial counts will be generated according to the distribution indicated by the \code{model} parameter.
#' @param SD The sequencing depth used to estimate capture efficiency; only used if the captureEfficient slot is not set in the \code{scaffoldParams} object.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simulateScaffold <- function(scaffoldParams, originalSCE, model = "p",
                             inputInitial=NULL, SD)
{
  if (is.null(inputInitial)) {
  initialCounts <- generateGeneCounts(numCells = scaffoldParams@numCells,
                                      mu = scaffoldParams@geneMeans,
                                      theta = scaffoldParams@geneTheta,
                                      type = model,
                                      degree = scaffoldParams@degree)
  rownames(initialCounts) <- scaffoldParams@genes
 } else {
  initialCounts = inputInitial
 }
  if (scaffoldParams@captureEfficiency[1] == -1) {
    print("about to estimate capture efficiency")
    capEfficiency <- estimateCaptureEff(Data = initialCounts,
                                       compareData = counts(originalSCE),
                                       SD = SD) # what is SD? seq depth??

  } else {
   capEfficiency <- scaffoldParams@captureEfficiency
  }

  print("finished capEfficiency")

  print("starting lysis")
  capturedMolecules <- lysisStep(round(initialCounts), capEfficiency)
  print("finished lysis")

  # Now we split between the different protocols

  if (scaffoldParams@protocol == "C1")
  {
    print("starting revt")
    efficiencyRT <- runif(scaffoldParams@numCells, .95, .95) # for now...
    capturedMolecules2 <- revtStep(capturedMolecules, efficiencyRT, Genes=scaffoldParams@genes)
    print("finished revt")

    print("starting preamplify")
    amplifiedMolecules <- preamplifyStep(capturedMolecules2,
                                         genes=scaffoldParams@genes,
                                         efficiencyPCR = scaffoldParams@firstAmpEfficiency,
                                         rounds = scaffoldParams@numFirstAmpCycles,
                                         typeAMP = scaffoldParams@typeOfAmp)
    print("finished preamplify")
    print("starting sequencing")
    finalCounts <- sequenceStepC1(amplifiedMolecules=amplifiedMolecules,
                                  pcntRange = scaffoldParams@percentRange,
                                  totalSD = scaffoldParams@totalSD,
                                  efficiencyPCR = scaffoldParams@secondAmpEfficiency,
                                  roundsPCR = scaffoldParams@numSecondAmpCycles,
                                  efficiencyTag = scaffoldParams@tagEfficiency)
    print("finished sequencing")
    SingleCellExperiment(assays = list(counts = finalCounts),
       metadata = list(initialSimCounts = initialCounts, capEfficiency = capEfficiency))
  }
  else if (scaffoldParams@protocol == "10x" || scaffoldParams@protocol == "10X")
  {
    print("starting sequencing")
    finalCounts <- sequenceStep10X(capturedMolecules, pcntRange = scaffoldParams@percentRange,
                                totalSD = scaffoldParams@totalSD,
                                efficiencyPCR = scaffoldParams@secondAmpEfficiency,
                                roundsPCR = scaffoldParams@numSecondAmpCycles,
                                genes=scaffoldParams@genes)
    print("finished sequencing")
    SingleCellExperiment(assays = list(counts = finalCounts))
  }

}
