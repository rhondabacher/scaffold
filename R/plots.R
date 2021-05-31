#' Generate summary plots of simulated data
#'
#' @param simulatedData An object of class SingleCellExperiment generated with the \code{simulateScaffold} function.
#' @param originalData The original SingleCellExperiment object used to initiate the simulation.
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom graphics lines par axis abline legend
#' @importFrom grDevices colorRampPalette
#' @export
makePlots <- function(simulatedData, originalData)
{
  simulatedCounts <- counts(simulatedData)
  originalCounts <- counts(originalData)
  Genes <- rownames(originalCounts)
  
  X <- data.frame( Depth = colSums(simulatedCounts), Species = "Simulated")
  Y <- data.frame( Depth = colSums(originalCounts), Species = "Original")
  longdata <- rbind(Y, X)
  par(mar=c(5,8,2,1), mgp = c(6, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
                    data = longdata,
                    xlab = "",
                    ylab = "Sequencing Depth", pal=c("cornflowerblue", "brown1"),
                    main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  
  X <- data.frame( Depth = colSums(simulatedCounts!=0) / nrow(simulatedCounts), Species = "Simulated")
  Y <- data.frame( Depth = colSums(originalCounts!=0) / nrow(simulatedCounts), Species = "Original")
  longdata <- rbind(Y, X)
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
                    data = longdata, ylim=c(0,1),
                    xlab = "",
                    ylab = "Cellular Detection Rate", pal=c("cornflowerblue", "brown1"),
                    main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  
  XX <- sample(Genes, 200)
  X1 <- apply(simulatedCounts[XX,], 1, function(x) sum(x!=0)) / dim(simulatedCounts)[2]
  X2 <- apply(originalCounts[XX,], 1, function(x) sum(x!=0)) / dim(originalCounts)[2]
  
  X <- data.frame( Depth =X1, Species = "Simulated")
  Y <- data.frame( Depth = X2, Species = "Original")
  longdata <- rbind(Y, X)
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
                    data = longdata,
                    xlab = "",
                    ylab = "Detection Rate (Gene)", pal=c("cornflowerblue", "brown1"),
                    main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  
  
  
  X1 <- log(apply(simulatedCounts[XX,], 1, function(x) mean(x))+1)
  X2 <- log(apply(originalCounts[XX,], 1, function(x) mean(x))+1)
  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
                    data = longdata,
                    xlab = "",
                    ylab = "log (mean+1)", pal=c("cornflowerblue", "brown1"),
                    main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  
  X1 <- log(apply(simulatedCounts[XX,], 1, function(x) sd(x))+1)
  X2 <- log(apply(originalCounts[XX,], 1, function(x) sd(x))+1)
  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
                    data = longdata,
                    xlab = "",
                    ylab = "log (sd+1)", pal=c("cornflowerblue", "brown1"),
                    main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  
  
  
}

