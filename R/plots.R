#' @importFrom SingleCellExperiment counts
#' @importFrom graphics lines par axis abline legend
#' @importFrom grDevices colorRampPalette
#' @importFrom SCnorm getSlopes
#' @importFrom stats median quantile sd ecdf
makePlots <- function(simulatedData, originalData)
{
  simulatedCounts <- counts(simulatedData)
  originalCounts <- counts(originalData)
  Genes <- rownames(originalCounts)

  par(mfrow=c(1,3), mar=c(5,5,2,1))
  options(mc.cores = 20)
  sim.slopes <- getSlopes(simulatedCounts, SeqDepth = colSums(simulatedCounts), Tau =.5, ditherCounts=FALSE)
  orig.slopes <- getSlopes(originalCounts, SeqDepth = colSums(originalCounts), Tau = .5, ditherCounts=FALSE)

  # Genes that have at least 10 non-zero in each dataset
  MedExpr <- log(apply(originalCounts[names(sim.slopes),], 1, function(x) median(x[x!=0])))
  splitby <- sort(MedExpr)
  grps <- length(splitby) / 10
  sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

  MedExpr <- log(apply(simulatedCounts[names(sim.slopes),], 1, function(x) median(x[x!=0])))
  splitby <- sort(MedExpr)
  grps <- length(splitby) / 10
  sreg_orig2 <- split(splitby, ceiling(seq_along(splitby) / grps))


  par(mfrow=c(2,2))
  makePlot(sim.slopes, sreg_orig2, "Simulated Data", ymax=4)
  makePlot(orig.slopes, sreg_orig1, "Original Data", ymax=1)
  dotPlot(sim.slopes, sreg_orig2, "Simulated Data")
  dotPlot(orig.slopes, sreg_orig1, "Original Data")

  X <- data.frame( Depth = colSums(simulatedCounts) / 1000000, Species = "Simulated")
  Y <- data.frame( Depth = colSums(originalCounts) / 1000000, Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "Sequencing Depth (millions)", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="Sequencing Depth (millions)",
       cex.axis=2, cex.lab=1.8, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X <- data.frame( Depth = colSums(simulatedCounts!=0) / nrow(simulatedCounts), Species = "Simulated")
  Y <- data.frame( Depth = colSums(originalCounts!=0) / nrow(simulatedCounts), Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata, ylim=c(0,1),
             xlab = "",
             ylab = "Cellular Detection Rate", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="Cellular Detection Rate",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(.2, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('topleft', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)

  XX <- sample(Genes, 100)
  X1 <- apply(simulatedCounts[XX,], 1, function(x) sum(x!=0)) / dim(simulatedCounts)[2]
  X2 <- apply(originalCounts[XX,], 1, function(x) sum(x!=0)) / dim(originalCounts)[2]

  X <- data.frame( Depth =X1, Species = "Simulated")
  Y <- data.frame( Depth = X2, Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "Dropout Rate", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="Dropout Rate",
       cex.axis=1.8, cex.lab=1.8, cex.main=1.8, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('topleft', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X1 <- apply(simulatedCounts, 1, function(x) sum(x!=0)) / dim(simulatedCounts)[2]
  X2 <- apply(originalCounts, 1, function(x) sum(x!=0)) / dim(originalCounts)[2]

  X <- data.frame( Depth =X1, Species = "Simulated")
  Y <- data.frame( Depth = X2, Species = "Original")
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="Dropout Rate",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('topleft', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)

  XX <- sample(Genes, 200)

  X1 <- apply(simulatedCounts[XX,], 1, function(x) sd(x))/ apply(simulatedCounts[XX,], 1, function(x) mean(x))
  X2 <- apply(originalCounts[XX,], 1, function(x) sd(x))/ apply(originalCounts[XX,], 1, function(x) mean(x))

  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "CV", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="CV",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X1 <- log(apply(simulatedCounts[XX,], 1, function(x) mean(x))+1)
  X2 <- log(apply(originalCounts[XX,], 1, function(x) mean(x))+1)
  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "log (mean+1)", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="log (mean+1)",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X1 <- log(apply(simulatedCounts[XX,], 1, function(x) sd(x))+1)
  X2 <- log(apply(originalCounts[XX,], 1, function(x) sd(x))+1)
  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "log (sd+1)", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="log (sd+1)",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X1 <- log(apply(simulatedCounts[XX,], 1, function(x) sd(x))+1)
  X2 <- log(apply(originalCounts[XX,], 1, function(x) sd(x[x < quantile(x,.99)]))+1)
  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "log (sd+1)", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="log (sd+1)",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X1 <- apply(simulatedCounts[XX,], 1, function(x) sd(x))/ apply(simulatedCounts[XX,], 1, function(x) mean(x))
  X2 <- apply(originalCounts[XX,], 1, function(x)
    sd(x[x < quantile(x,.95)]))/
    apply(originalCounts[XX,], 1, function(x) mean(x[x < quantile(x,.95)]))

  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "CV", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=2, cex.axis=2,cex.names=2)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="CV",
       cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)

  X1 <- apply(simulatedCounts[XX,], 1, function(x) sd(x))
  X2 <- apply(originalCounts[XX,], 1, function(x) sd(x[x < quantile(x,.95)]))
  useg <- names(which(X1<Inf & X2 < Inf))
  X <- data.frame( Depth =X1[useg], Species = "Simulated")
  Y <- data.frame( Depth = X2[useg], Species = "Original")
  longdata <- rbind(Y, X)
  par(mfrow=c(1,2))
  par(mar=c(5,7,2,1), mgp = c(4, .5, 0))
  yarrr::pirateplot(formula = Depth ~ Species,
             data = longdata,
             xlab = "",
             ylab = "CV", pal=c("cornflowerblue", "brown1"),
             main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.6, cex.axis=1.5,cex.names=1.6)
  par(mar=c(5,7,2,1), mgp = c(3, 1, 0))
  plot(ecdf(X$Depth), col="brown1", main="", xlab="CV",
       cex.axis=1.5, cex.lab=1.6, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
  plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
  legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)
}




makePlot<- function(SLOPES, sreg, TYPE, ymax=1) {
  colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = length(sreg))
  plot(density(na.omit(SLOPES), from=min(SLOPES, na.rm=T), to=max(SLOPES, na.rm=T), adjust=1),
       xlab="Slope", ylab="Density",  main=TYPE,
       cex.lab=1.5, cex.main=1.5, cex.axis=1.6,xlim=c(-2,3),
       ylim=c(0,ymax), lwd=3, col="white", xaxt='n', yaxt='n', bty='n')
  for (i in rev(which(sapply(sreg, length) > 3))) {
    topgenes<-names(sreg[[i]])
    lines(density(na.omit(SLOPES[topgenes]), from=-3, to=3, adjust=1), lwd=3, col=colors[i])
  }
  abline(v=0, lwd=2, col="black")
  axis(side = 1, at=seq(-2,3, by=1), lwd.ticks=2, lwd = 2, cex.axis=1.2,font=1)
  axis(side = 2, at=seq(0,ymax,by=.5), lwd = 2, lwd.ticks=2, cex.axis=1.2,font=1)
  abline(v=1, lwd=2, col="black")

}

makeCDP <- function(counts, NAME) {

  gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

  MedExpr <- log(apply(counts[names(gslopes),], 1, function(x) median(x[x!=0])))
  splitby <- sort(MedExpr)
  grps <- length(splitby) / 10
  sreg <- split(splitby, ceiling(seq_along(splitby) / grps))
  makePlot(gslopes, NAME, sreg)
}

makeCDPdot <- function(counts, NAME) {

  gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

  MedExpr <- log(apply(counts[names(gslopes),], 1, function(x) median(x[x!=0])))
  splitby <- sort(MedExpr)
  grps <- length(splitby) / 10
  sreg <- split(splitby, ceiling(seq_along(splitby) / grps))
  dotPlot(gslopes,sreg, NAME)
}


#' @importFrom stats na.omit density
dotPlot <- function(SLOPES, sreg, NAME) {
  colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = length(sreg))
  Mode <- c()
  DensH <- c()
  for (i in which(sapply(sreg, length) > 3)) {
    useg <- names(sreg[[i]])
    rqdens <- density(na.omit(SLOPES[useg]), from=-3, to=3)
    peak <- which.max(rqdens$y)
    Mode[i] <- rqdens$x[peak]
    DensH[i] <- rqdens$y[peak]
  }

  par(mar=c(5,5,2,1))
  plot(Mode, 1:length(sreg), xlim=c(-.5,1.5), pch=19, lwd=3, col=colors, xlab="Slope", bty='n',xaxt='n', yaxt='n',
       cex=1.5, cex.lab=1.8, cex.axis=1.5, cex.main=1.5, ylab="Expression Group", main=NAME)
  axis(side = 2, at=1:length(sreg), lwd.ticks=2, lwd = 2, cex.axis=1.2,font=1)
  axis(side = 1, lwd = 2, lwd.ticks=2, cex.axis=1.2,font=1)
  abline(v=1, lty=2, lwd=3)

  # Modestat <- sum((Mode[3:10] - 1)^2)
  # return(Modestat)

}
