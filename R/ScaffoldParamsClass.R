setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))

#' S4 object for scaffold simulation parameters
#'
#' @export
setClass("ScaffoldParams",
         representation(numCells = "numeric",
                        numGenes = "numeric",
                        geneMeans = "numeric",
                        geneTheta = "numericOrNULL",
												totalTranscripts = "numeric",
                        genes = "character",
                        captureEfficiency = "numeric",
												efficiencyRT = "numericOrNULL",
                        typeOfAmp = "character",
                        numFirstAmpCycles = "numeric",
                        numSecondAmpCycles = "numeric",
                        firstAmpEfficiency = "numeric",
                        secondAmpEfficiency = "numeric",
                        tagEfficiency = "numeric",
                        degree = "numeric",
                        percentRange = "numeric",
                        protocol = "character",
                        totalSD = "numeric",
												useUMI = "logical",
												model = "character",
												usePops = "listOrNULL",
												useDynamic = "logical",
												propDynamic = "numeric"))

#' @importFrom utils head
setMethod("show", signature(object = "ScaffoldParams"),
          function(object){
            cat("protocol:", object@protocol, "\n")
            cat("numCells:", object@numCells,"\n")
            cat("numGenes:",object@numGenes, "\n")
						cat("totalTranscripts:",object@totalTranscripts, "\n")
            cat("geneMeans:", head(object@geneMeans), "...\n")
            cat("geneTheta:", head(object@geneTheta), "...\n")
            cat("genes:", head(object@genes), "...\n")
            cat("captureEfficiency:", head(object@captureEfficiency), "\n")
						cat("efficiencyRT:", head(object@efficiencyRT), "\n")
            cat("typeOfAmp:", object@typeOfAmp, "\n")
            cat("numFirstAmpCycles:", object@numFirstAmpCycles, "\n")
            cat("numSecondAmpCycles:", object@numSecondAmpCycles, "\n")
            cat("firstAmpEfficiency:", head(object@firstAmpEfficiency), "...\n")
            cat("secondAmpEfficiency:", head(object@secondAmpEfficiency), "...\n")
            cat("tagEfficiency:", head(object@tagEfficiency), "...\n")
            cat("degree:", object@degree, "\n")
            cat("percentRange:", object@percentRange, "\n")
            cat("totalSD:", object@totalSD, "\n")
						cat("useUMI:", object@useUMI, "\n")
						cat("model:", object@model, "\n")
						print(paste0("usePops:", head(object@usePops), "\n"))

          })
