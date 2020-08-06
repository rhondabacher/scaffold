#' Run a scaffold simulation.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simulateScaffold <- function(scaffoldParams, originalSCE, model = "p", inputInitial=NULL, SD)
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
