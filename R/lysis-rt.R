                                        # this is step 2
#' @importFrom wrswoR sample_int_rej
lysisStep <- function(Data, efficiencyCell = NULL, efficiencyGene = NULL){

  if (is.null(efficiencyGene) == TRUE) {
    efficiencyGene <- rep(1, nrow(Data))
  }

  efficiencyGene <- efficiencyGene / sum(efficiencyGene)
  names(efficiencyGene) <- rownames(Data)

  Genes <- rownames(Data)
  allMolecules <- apply(Data, 2, function(Cell) rep(Genes, Cell))# list

  capturedMolecules <- lapply(1:length(allMolecules), function(x) {
  countAll <- Data[,x]
  efficiencyG <- efficiencyGene[names(countAll)]
  geneProbs <- rep(efficiencyG, countAll)
  geneProbs <- geneProbs / sum(geneProbs)

  totalM <- abs(round(efficiencyCell[x]*length(allMolecules[[x]])))

  sampledM <- sample_int_rej(n = length(allMolecules[[x]]), size = totalM, prob = geneProbs)
  sampledM <- allMolecules[[x]][sampledM]

  return(sampledM)
  })

}

revtStep <- function(Data, efficiencyCell = NULL, efficiencyGene = NULL,
                     Genes=NULL){

  Size <- unlist(lapply(Data, function(x) sum(table(x))))



  if(is.null(efficiencyGene) == TRUE) {
    efficiencyGene <- rep(1, length(Genes))
  }
  efficiencyGene <- efficiencyGene / sum(efficiencyGene)
  names(efficiencyGene) <- Genes

  capturedMolecules <- lapply(1:length(Data), function(x) {
    countAll <- table(Data[[x]])

    efficiencyG <- efficiencyGene[names(countAll)]
    geneProbs <- rep(efficiencyG, countAll)
    geneProbs <- geneProbs / sum(geneProbs)

    totalM <- abs(round(efficiencyCell[x]*length(Data[[x]])))
    sampledM <- sample_int_rej(n = length(Data[[x]]), size = totalM, prob = geneProbs)
    sampledM <- Data[[x]][sampledM]

    # print(x)
    return(sampledM)
  })
}
