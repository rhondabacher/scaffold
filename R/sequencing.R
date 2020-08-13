sequenceStepC1 <- function(amplifiedMolecules, pcntRange=0, totalSD=50000000,
                         efficiencyPCR, roundsPCR, efficiencyTag=NULL) {

  amplifiedMoleculesQuant <- quantCells(amplifiedMolecules, pcntRange=pcntRange)

  Size <- unlist(lapply(amplifiedMoleculesQuant, sum))

  taggedMolecules <- lapply(1:length(amplifiedMoleculesQuant), function(x) {
    geneProbs_all <- log(amplifiedMoleculesQuant[[x]] / sum(as.numeric(amplifiedMoleculesQuant[[x]])))
    totalT <- round(efficiencyTag[x]*sum(as.numeric(amplifiedMoleculesQuant[[x]])))
    tagdM <- try(rmultinom(n=1, size=totalT, prob=exp(geneProbs_all)), silent=T)
    countsT <- as.vector(tagdM)
    names(countsT) <- rownames(tagdM)
    return(countsT)
  })

  taggedMolecules <- lapply(taggedMolecules, function(x) x * 10)

  quantifiedMolecules <- amplifyStep(taggedMolecules, genes=names(taggedMolecules[[1]]), efficiencyPCR, roundsPCR, protocol = "C1")

  amplifiedMoleculesQuant_all_list <- lapply(1:length(quantifiedMolecules), function(x) {
    y = quantifiedMolecules[[x]
    names(y) <- paste0(names(quantifiedMolecules[[x]]),"_Cell_", x)
    return(y)})

  amplifiedMoleculesQuant_all <- do.call(c, amplifiedMoleculesQuant_all_list)
  geneProbs_all <- log(amplifiedMoleculesQuant_all / sum(as.numeric(amplifiedMoleculesQuant_all)))

  counts <- try(rmultinom(n=1, size=totalSD, prob=exp(geneProbs_all)), silent=T)

  cell_split <- matrix(counts, nrow=length(amplifiedMoleculesQuant[[1]]), ncol=length(amplifiedMoleculesQuant), byrow=FALSE)
  rownames(cell_split) <- names(quantifiedMolecules[[1]])
  colnames(cell_split) <- paste0("Cell_", 1:length(amplifiedMoleculesQuant))

  return(cell_split)

}

# under construction
sequenceStep10X <- function(amplifiedMolecules, pcntRange=0, totalSD=50000000,
                            efficiencyPCR, roundsPCR, genes)
{
  print("rearranging 10x data")

  ampMol.updated <- lapply(1:length(amplifiedMolecules), function(x) paste0(amplifiedMolecules[[x]], "_Cell", x))

  longVectorCounts <- table(unlist(ampMol.updated))

  quantifiedMolecules <- amplifyStep(longVectorCounts, genes=genes, efficiencyPCR, roundsPCR, protocol = "10X")


  amplifiedMoleculesQuant_all <- do.call(c, quantifiedMolecules)


  geneProbs_all <- log(amplifiedMoleculesQuant_all / sum(as.numeric(amplifiedMoleculesQuant_all)))


  counts <- try(rmultinom(n=1, size=totalSD, prob=exp(geneProbs_all)), silent=T)


  cell_split <- matrix(counts, nrow=length(quantifiedMolecules[[1]]), ncol=length(quantifiedMolecules), byrow=FALSE)
  rownames(cell_split) <- genes
  colnames(cell_split) <- paste0("Cell_", 1:length(quantifiedMolecules))

  return(cell_split)
}

