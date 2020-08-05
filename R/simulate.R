#' @importFrom SingleCellExperiment SingleCellExperiment
simulateSplash <- function(splashParams, originalSCE, model = "p", inputInitial=NULL, SD)
{
  if (is.null(inputInitial)) {
  initialCounts <- generateGeneCounts(numCells = splashParams@numCells,
                                      mu = splashParams@geneMeans,
                                      theta = splashParams@geneTheta,
                                      type = model,
                                      degree = splashParams@degree)
  rownames(initialCounts) <- splashParams@genes
 } else {
  initialCounts = inputInitial
 }
  if (splashParams@captureEfficiency[1] == -1) {
    print("about to estimate capture efficiency")
    capEfficiency <- estimateCaptureEff(Data = initialCounts,
                                       compareData = counts(originalSCE),
                                       SD = SD) # what is SD? seq depth??

  } else {
   capEfficiency <- splashParams@captureEfficiency
  }

  print("finished capEfficiency")

  print("starting lysis")
  capturedMolecules <- lysisStep(round(initialCounts), capEfficiency)
  print("finished lysis")

  # Now we split between the different protocols

  if (splashParams@protocol == "C1")
  {
    print("starting revt")
    efficiencyRT <- runif(splashParams@numCells, .95, .95) # for now...
    capturedMolecules2 <- revtStep(capturedMolecules, efficiencyRT, Genes=splashParams@genes)
    print("finished revt")

    print("starting preamplify")
    amplifiedMolecules <- preamplifyStep(capturedMolecules2,
                                         genes=splashParams@genes,
                                         efficiencyPCR = splashParams@firstAmpEfficiency,
                                         rounds = splashParams@numFirstAmpCycles,
                                         typeAMP = splashParams@typeOfAmp)
    print("finished preamplify")
    print("starting sequencing")
    finalCounts <- sequenceStepC1(amplifiedMolecules=amplifiedMolecules,
                                  pcntRange = splashParams@percentRange,
                                  totalSD = splashParams@totalSD,
                                  efficiencyPCR = splashParams@secondAmpEfficiency,
                                  roundsPCR = splashParams@numSecondAmpCycles,
                                  efficiencyTag = splashParams@tagEfficiency)
    print("finished sequencing")
    SingleCellExperiment(assays = list(counts = finalCounts),
       metadata = list(initialSimCounts = initialCounts, capEfficiency = capEfficiency))
  }
  else if (splashParams@protocol == "10x" || splashParams@protocol == "10X")
  {
    print("starting sequencing")
    finalCounts <- sequenceStep10X(capturedMolecules, pcntRange = splashParams@percentRange,
                                totalSD = splashParams@totalSD,
                                efficiencyPCR = splashParams@secondAmpEfficiency,
                                roundsPCR = splashParams@numSecondAmpCycles,
                                genes=splashParams@genes)
    print("finished sequencing")
    SingleCellExperiment(assays = list(counts = finalCounts))
  }

}
