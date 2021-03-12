#' Run a scaffold simulation.
#'
#' @param scaffoldParams An object of class ScaffoldParams.
#' @param originalSCE The SingleCellExperiment used to create the scaffoldParams object.
#' @param model The probability distribution used to generate the initial gene counts. Defaults to 'p' for the Poisson distribution, but can also be set to 'nb' to use the Negative Binomial distribution.
#' @param inputInitial A optional matrix of initial gene counts which should have the same dimension indicated in the \code{scaffoldParams} parameter. If left NULL, the initial counts will be generated according to the distribution indicated by the \code{model} parameter.
#' @param SD The standard deviation used to estimate capture efficiency; only used if the captureEfficiency slot is not set in the \code{scaffoldParams} object.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simulateScaffold <- function(scaffoldParams, originalSCE, model = "p",
                             inputInitial=NULL, SD = 0.02)
{
	
	numCells <- sum(scaffoldParams@numCells)
	cellPopulation <- rep(1:length(scaffoldParams@numCells), scaffoldParams@numCells)
	
	if (!is.null(scaffoldParams@sepPops[[1]])) {
	  initialCounts <- generateGeneCounts(numCells = numCells,
	                                      mu = scaffoldParams@geneMeans,
	                                      theta = scaffoldParams@geneTheta,
	                                      type = model,
	                                      degree = scaffoldParams@degree)
	  rownames(initialCounts) <- scaffoldParams@genes
    
	  cellSplit <- split(1:numCells, f=cellPopulation)
  
	  temp <- lapply(2:length(cellSplit), function(x) {
    
	    pop_temp <- initialCounts[,cellSplit[[x]]]
	    numSamp <- nrow(initialCounts) * scaffoldParams@sepPops$propGenes[(x-1)]
	    selectGenes <- sample(1:nrow(initialCounts), numSamp)			
	    fc_genes <- rnorm(length(selectGenes), mean = scaffoldParams@sepPops$fc[(x-1)], sd=.4) ## Get these from somewhere!!
	    fc_mat <- sapply(fc_genes, function(y) {
	      rnorm(ncol(pop_temp), mean = y, sd = .1) #so not every value multipled exactly the same
	    })
	    pop_temp[selectGenes,] <- pop_temp[selectGenes,] * t(fc_mat)
    
    
	    initialCounts[,cellSplit[[x]]] <- pop_temp
	    return()
	  })
  
	} else if(is.null(scaffoldParams@sepPops[[1]])) {
	  if (is.null(inputInitial)) {
	    initialCounts <- generateGeneCounts(numCells = numCells,
	                                        mu = scaffoldParams@geneMeans,
	                                        theta = scaffoldParams@geneTheta,
	                                        type = model,
	                                        degree = scaffoldParams@degree)
	    rownames(initialCounts) <- scaffoldParams@genes
	  } else {
	    initialCounts = inputInitial
	  }
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

  print("starting lysis and rev tran")
  efficiencyRT <- runif(numCells, .95, .95) # for now...
  capturedMolecules <- captureStep(round(initialCounts), captureEffCell=capEfficiency, 
                                   rtEffCell = efficiencyRT, useUMI = scaffoldParams@useUMI)
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
       metadata = list(initialSimCounts = initialCounts, capEfficiency = capEfficiency, cellPopulation = cellPopulation))
  }
  else if (scaffoldParams@protocol == "10x" || scaffoldParams@protocol == "10X")
  {
    print("starting sequencing")
    finalCounts <- sequenceStep10X(capturedMolecules, 
                                totalSD = scaffoldParams@totalSD,
                                efficiencyPCR = scaffoldParams@secondAmpEfficiency,
                                roundsPCR = scaffoldParams@numSecondAmpCycles,
                                genes = scaffoldParams@genes,
																useUMI = scaffoldParams@useUMI)
    print("finished sequencing")
    SingleCellExperiment(assays = list(counts = finalCounts$counts, umi_counts=finalCounts$umi_counts),
       metadata = list(initialSimCounts = initialCounts, capEfficiency = capEfficiency, cellPopulation = cellPopulation))
  }

}
