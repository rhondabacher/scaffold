#' Run a Scaffold simulation.
#'
#' @param scaffoldParams An object of class ScaffoldParams. Generated using the Scaffold::estimateScaffoldParameters function.
#' @param originalSCE The SingleCellExperiment used to create the scaffoldParams object.
#' @param inputInitial A optional matrix of initial gene counts which should have the same dimension indicated in the \code{scaffoldParams} parameter. If left NULL, the initial counts will be generated according to the distribution indicated by the \code{model} parameter. This is mainly used in simulations to isolate the effects of each step without regenerating a new initial mRNA counts.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simulateScaffold <- function(scaffoldParams, originalSCE, inputInitial=NULL)
{
  
  numCells <- sum(scaffoldParams@numCells)
  cellPopulation <- rep(1:length(scaffoldParams@numCells), scaffoldParams@numCells)
  
  geneStatus <- rep(NA, scaffoldParams@numGenes)
  # Simulating dynamic populations
  if (!is.null(scaffoldParams@useDynamic[[1]])) {
    dynamicsim <- generateDynamicGeneCounts(numCells = numCells, 
                                            mu = scaffoldParams@geneMeans, 
                                            dynamicParams = scaffoldParams@useDynamic)
    initialCounts <- dynamicsim[[1]]
    geneStatus <- dynamicsim[[2]]
  } else {
    if (!is.null(scaffoldParams@usePops[[1]])) {
      cellSplit <- split(1:numCells, f=cellPopulation)
      
      allCounts <- lapply(1:length(cellSplit), function(x){
        
        means <- scaffoldParams@geneMeans
        numSamp <- length(means) * scaffoldParams@usePops$propGenes[x]
        if (numSamp > 0) {
          selectGenes <- sample(1:length(means), numSamp)
          fc_genes <- abs(rnorm(length(selectGenes), mean = scaffoldParams@usePops$fc_mean[x], sd=scaffoldParams@usePops$fc_sd[x]))
          flipfc <- sample(1:length(fc_genes), length(fc_genes) / 2)
          fc_genes[flipfc] <- 1/ fc_genes[flipfc]
          
          means[selectGenes] <- means[selectGenes] * fc_genes
        }
        generateCnts <- generateGeneCounts(numCells = scaffoldParams@numCells[x],
                                           mu = means,
                                           popHet = scaffoldParams@popHet)
        rownames(generateCnts) <- scaffoldParams@genes
        return(generateCnts)
      })
      initialCounts <- do.call(cbind, allCounts)
      
    } else if (is.null(scaffoldParams@usePops[[1]])) {
      if (is.null(inputInitial)) {
        initialCounts <- generateGeneCounts(numCells = numCells,
                                            mu = scaffoldParams@geneMeans,
                                            popHet = scaffoldParams@popHet)
        rownames(initialCounts) <- scaffoldParams@genes
      } else {
        initialCounts = inputInitial
      }
    }
  }
  
  if (is.null(scaffoldParams@captureEfficiency)) {
    print("Estimating capture efficiency...")
    capEfficiency <- estimateCaptureEff(Data = initialCounts,
                                        compareData = counts(originalSCE),
                                        protocol = scaffoldParams@protocol,
                                        fromUMI = scaffoldParams@sceUMI)
  } else {
    capEfficiency <- scaffoldParams@captureEfficiency
  }
  print("Finished estimating capture efficiency!")
  
  print("Starting capture step (lysis and reverse transcription)...")
  capturedMolecules <- captureStep(round(initialCounts), 
                                   captureEffCell = capEfficiency, 
                                   rtEffCell = scaffoldParams@efficiencyRT, 
                                   useUMI = scaffoldParams@useUMI)
  print("Finished capture step!")
  
  # Now we split between the different protocols
  
  if (scaffoldParams@protocol == "FULLLENGTH") {
    print("Starting preamplify step...")
    amplifiedMolecules <- preamplifyStep(capturedMolecules = capturedMolecules,
                                         genes = scaffoldParams@genes,
                                         efficiencyPCR = scaffoldParams@preAmpEfficiency,
                                         rounds = scaffoldParams@numPreAmpCycles,
                                         typeAMP = scaffoldParams@typeOfAmp,
                                         useUMI = scaffoldParams@useUMI)
    print("Finished preamplify!")
    print("Starting library prep and sequencing...")
    finalCounts <- sequenceStepC1(amplifiedMolecules=amplifiedMolecules,
                                  equalizationAmount = scaffoldParams@equalizationAmount,
                                  totalDepth = scaffoldParams@totalDepth,
                                  efficiencyPCR = scaffoldParams@ampEfficiency,
                                  roundsPCR = scaffoldParams@numAmpCycles,
                                  efficiencyTag = scaffoldParams@tagEfficiency,
                                  genes = scaffoldParams@genes,
                                  useUMI = scaffoldParams@useUMI)
    print("Finished sequencing and data formatting!")
    
    if(is.null(finalCounts$umi_counts)) {
      finalCounts$umi_counts <- matrix(NA, nrow=nrow(finalCounts$counts), ncol=ncol(finalCounts$counts))
    }
    
    returnsce <- SingleCellExperiment(assays = list(counts = finalCounts$counts, umi_counts=finalCounts$umi_counts),
                         metadata = list(initialSimCounts = initialCounts),
                         colData = data.frame(capEfficiency = capEfficiency, cellPopulation = factor(cellPopulation)),
                         rowData = data.frame(geneStatus = geneStatus))
  }
  
  if (scaffoldParams@protocol == "DROPLET") {
    print("Starting library prep and sequencing...")
    finalCounts <- sequenceStep10X(capturedMolecules, 
                                   totalDepth = scaffoldParams@totalDepth,
                                   efficiencyPCR = scaffoldParams@ampEfficiency,
                                   roundsPCR = scaffoldParams@numAmpCycles,
                                   genes = scaffoldParams@genes,
                                   efficiencyTag = scaffoldParams@tagEfficiency,
                                   useUMI = scaffoldParams@useUMI)
    print("Finished sequencing and data formatting!")
    returnsce <- SingleCellExperiment(assays = list(counts = finalCounts$counts, umi_counts=finalCounts$umi_counts),
                         metadata = list(initialSimCounts = initialCounts),
                         colData = data.frame(capEfficiency = capEfficiency, cellPopulation = factor(cellPopulation)),
                         rowData = data.frame(geneStatus = geneStatus))
  }
  return(returnsce)
}
