setClassUnion("numericOrNULL", members = c("numeric", "NULL"))
setClassUnion("listOrNULL", members = c("list", "NULL"))

#' S4 object for Scaffold simulation parameters
#'
#' @export
setClass("ScaffoldParams",
         representation(sceUMI = "logical",
					 							numCells = "numeric",
                        numGenes = "numeric",
                        geneMeans = "numeric",
												totalTranscripts = "numeric",
                        genes = "character",
												protocol = "character",
												useUMI = "logical",
												popHet = "numeric",
												geneEfficiency = "numericOrNULL",
                        captureEfficiency = "numericOrNULL",
												efficiencyRT = "numericOrNULL",
                        typeOfAmp = "character",
												numPreAmpCycles = "numeric",
                        numAmpCycles = "numeric",
                        preAmpEfficiency = "numeric",
                        ampEfficiency = "numeric",
                        tagEfficiency = "numericOrNULL",
                        equalizationAmount = "numeric",
                        totalDepth = "numeric",
                        usePops = "listOrNULL",
												useDynamic = "listOrNULL"))

#' @importFrom utils head
setMethod("show", signature(object = "ScaffoldParams"),
          function(object){
            cat("protocol:", object@protocol, "\n")
            cat("numCells:", object@numCells,"\n")
            cat("numGenes:",object@numGenes, "\n")
						cat("totalTranscripts:",object@totalTranscripts, "\n")
            cat("geneMeans:", head(object@geneMeans), "...\n")
            cat("genes:", head(object@genes), "...\n")
						cat("popHet:", object@popHet, "\n")
						cat("useUMI:", object@useUMI, "\n")
            cat("captureEfficiency:", head(object@captureEfficiency), "\n")
						cat("efficiencyRT:", head(object@efficiencyRT), "\n")
            cat("typeOfAmp:", object@typeOfAmp, "\n")
            cat("numPreAmpCycles:", object@numFirstAmpCycles, "\n")
            cat("numAmpCycles:", object@numSecondAmpCycles, "\n")
            cat("preAmpEfficiency:", head(object@firstAmpEfficiency), "...\n")
            cat("ampEfficiency:", head(object@secondAmpEfficiency), "...\n")
            cat("tagEfficiency:", head(object@tagEfficiency), "...\n")
            cat("equalizationAmount:", object@equalizationAmount, "\n") 
            cat("totalDepth:", object@totalDepth, "\n")
            
						print(paste0("usePops:", head(object@usePops), "\n"))
						print(paste0("useDynamic:", head(object@useDynamic), "\n"))
          })
