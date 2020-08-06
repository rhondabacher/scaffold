setClass("ScaffoldParams",
         representation(numCells = "numeric",
                        numGenes = "numeric",
                        geneMeans = "numeric",
                        geneTheta = "numeric",
                        genes = "character",
                        captureEfficiency = "numeric",
                        typeOfAmp = "character",
                        numFirstAmpCycles = "numeric",
                        numSecondAmpCycles = "numeric",
                        firstAmpEfficiency = "numeric",
                        secondAmpEfficiency = "numeric",
                        tagEfficiency = "numeric",
                        degree = "numeric", # should we keep this?
                        percentRange = "numeric",
                        protocol = "character",
                        totalSD = "numeric"))

#' @importFrom utils head
setMethod("show", signature(object = "ScaffoldParams"),
          function(object){
            cat("protocol:", object@protocol, "\n")
            cat("numCells:", object@numCells,"\n")
            cat("numGenes:",object@numGenes, "\n")
            cat("geneMeans:", head(object@geneMeans), "...\n")
            cat("geneTheta:", head(object@geneTheta), "...\n")
            cat("genes:", head(object@genes), "...\n")
            cat("captureEfficiency:", head(object@captureEfficiency), "\n")
            cat("typeOfAmp:", object@typeOfAmp, "\n")
            cat("numFirstAmpCycles:", object@numFirstAmpCycles, "\n")
            cat("numSecondAmpCycles:", object@numSecondAmpCycles, "\n")
            cat("firstAmpEfficiency:", head(object@firstAmpEfficiency), "...\n")
            cat("secondAmpEfficiency:", head(object@secondAmpEfficiency), "...\n")
            cat("tagEfficiency:", head(object@tagEfficiency), "...\n")
            cat("degree:", object@degree, "\n")
            cat("percentRange:", object@percentRange, "\n")
            cat("totalSD:", object@totalSD, "\n")

          })
