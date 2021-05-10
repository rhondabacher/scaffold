#' Run a scaffold simulation.
#'
#' @param scaffoldParams An object of class ScaffoldParams.
#' @param originalSCE The SingleCellExperiment used to create the scaffoldParams object.
#' @param inputInitial A optional matrix of initial gene counts which should have the same dimension indicated in the \code{scaffoldParams} parameter. If left NULL, the initial counts will be generated according to the distribution indicated by the \code{model} parameter.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simulateScaffold <- function(scaffoldParams, originalSCE, inputInitial=NULL)
{
	
	numCells <- sum(scaffoldParams@numCells)
	cellPopulation <- rep(1:length(scaffoldParams@numCells), scaffoldParams@numCells)
	
	if (scaffoldParams@useDynamic == TRUE) {
		initialCounts <- generateDynamicGeneCounts(numCells = numCells, 
																							mu = scaffoldParams@geneMeans, 
																							propDynamic = scaffoldParams@propDynamic)
	rownames(initialCounts) <- scaffoldParams@genes
	} else {
	if (!is.null(scaffoldParams@usePops[[1]])) {
		cellSplit <- split(1:numCells, f=cellPopulation)
    
		    allCounts <- lapply(1:length(cellSplit), function(x){
      
		      means <- scaffoldParams@geneMeans
		      numSamp <- length(means) * scaffoldParams@usePops$propGenes[(x)]
		      if (numSamp > 0) {
		        selectGenes <- sample(1:length(means), numSamp)
		        fc_genes <- abs(rnorm(length(selectGenes), mean = scaffoldParams@usePops$fc[(x)], sd=.4))
		        flipfc <- sample(1:length(fc_genes), length(fc_genes) / 2)
		        fc_genes[flipfc] <- 1/ fc_genes[flipfc]
      
		        means[selectGenes] <- means[selectGenes] * fc_genes
		      }
		      generateCnts <- generateGeneCounts(numCells = scaffoldParams@numCells[x],
		                                          mu = means,
		                                          theta = scaffoldParams@geneTheta,
		                                          type = scaffoldParams@model,
		                                          degree = scaffoldParams@degree)
		      rownames(generateCnts) <- scaffoldParams@genes
		      return(generateCnts)
		    })
    
		    initialCounts <- do.call(cbind, allCounts)

  
	} else if(is.null(scaffoldParams@usePops[[1]])) {
	  if (is.null(inputInitial)) {
	    initialCounts <- generateGeneCounts(numCells = numCells,
	                                        mu = scaffoldParams@geneMeans,
	                                        theta = scaffoldParams@geneTheta,
	                                        type = scaffoldParams@model,
	                                        degree = scaffoldParams@degree)
	    rownames(initialCounts) <- scaffoldParams@genes
	  } else {
	    initialCounts = inputInitial
	  }
	}
}
	
  

 
 
  if (scaffoldParams@captureEfficiency[1] == -1) {
    print("about to estimate capture efficiency")
    capEfficiency <- estimateCaptureEff(Data = initialCounts,
                                       compareData = counts(originalSCE),
                                       protocol = scaffoldParams@protocol,
                                       fromUMI = scaffoldParams@fromUMI)

  } else {
   capEfficiency <- scaffoldParams@captureEfficiency
  }

  print("finished capEfficiency")

  print("starting lysis and rev tran")
  capturedMolecules <- captureStep(round(initialCounts), 
																		captureEffCell=capEfficiency, 
                                    rtEffCell = scaffoldParams@efficiencyRT, 
																		useUMI = scaffoldParams@useUMI)
  print("finished lysis and revt")

  # Now we split between the different protocols

  if (scaffoldParams@protocol == "C1")
  {
    print("starting preamplify")
    amplifiedMolecules <- preamplifyStep(capturedMolecules,
                                         genes=scaffoldParams@genes,
                                         efficiencyPCR = scaffoldParams@firstAmpEfficiency,
                                         rounds = scaffoldParams@numFirstAmpCycles,
                                         typeAMP = scaffoldParams@typeOfAmp,
                                         useUMI = scaffoldParams@useUMI)
    print("finished preamplify")
    print("starting sequencing")
    finalCounts <- sequenceStepC1(amplifiedMolecules=amplifiedMolecules,
                                  pcntRange = scaffoldParams@percentRange,
                                  totalSD = scaffoldParams@totalSD,
                                  efficiencyPCR = scaffoldParams@secondAmpEfficiency,
                                  roundsPCR = scaffoldParams@numSecondAmpCycles,
                                  efficiencyTag = scaffoldParams@tagEfficiency,
                                  genes = scaffoldParams@genes,
                                  useUMI = scaffoldParams@useUMI)
    print("finished sequencing")
		if(is.null(finalCounts$umi_counts)) {
			finalCounts$umi_counts <- matrix(NA, nrow=nrow(finalCounts$counts), ncol=ncol(finalCounts$counts))
		}
			
    SingleCellExperiment(assays = list(counts = finalCounts$counts, umi_counts=finalCounts$umi_counts),
       metadata = list(initialSimCounts = initialCounts),
			 colData = data.frame(capEfficiency = capEfficiency, cellPopulation = factor(cellPopulation)))
  }
  else if (scaffoldParams@protocol == "10x" || scaffoldParams@protocol == "10X")
  {
    print("starting sequencing")
    finalCounts <- sequenceStep10X(capturedMolecules, 
                                totalSD = scaffoldParams@totalSD,
                                efficiencyPCR = scaffoldParams@secondAmpEfficiency,
                                roundsPCR = scaffoldParams@numSecondAmpCycles,
                                genes = scaffoldParams@genes,
																efficiencyTag = scaffoldParams@tagEfficiency,
																useUMI = scaffoldParams@useUMI)
    print("finished sequencing")
    SingleCellExperiment(assays = list(counts = finalCounts$counts, umi_counts=finalCounts$umi_counts),
       metadata = list(initialSimCounts = initialCounts),
			 colData = data.frame(capEfficiency = capEfficiency, cellPopulation = factor(cellPopulation)))
  }

}
