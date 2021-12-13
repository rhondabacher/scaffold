#' Generate initial counts
#'
#' Generates the transcripts in the cell. 
#'
#' @param numCells The number of cells to generate counts for.
#' @param mu a vector containing the mean of each gene to simulate counts for.
#' @param popHet the degree of heterogeneity within the cell population.
#'	Should be a vector of length two specifying the bounds.
#'
#' @return Output is matrix of counts, genes on the rows and cells
#'   on the columns.
#'
#' @importFrom stats rpois
#' @importFrom Rfast Sort
generateGeneCounts <- function(numCells, mu, popHet) {
	
    R <- matrix(sapply(1:length(mu), function(x) rpois(numCells, mu[x])), nrow=numCells)

    R <- t(R*Rfast::Sort(runif(nrow(R), popHet[1], popHet[2])))
  
  return(R)
}

#' Generate a dynamic cell population
## The motivation for this code is from the simstudy package.
#' @inheritParams generateGeneCounts
#' @inheritParams estimateScaffoldParameters
#' @import splines
generateDynamicGeneCounts <- function(numCells, mu, dynamicParams) {
  

    selectGenes <- sample(names(mu), ceiling(dynamicParams$propGenes*length(mu)))
    otherGenes <- setdiff(names(mu), selectGenes)

  R0 <- matrix(sapply(otherGenes, function(x) rpois(numCells, mu[x])), nrow=numCells)

  if (is.null(dynamicParams$degree)) dynamicParams$degree <- 2

  ngenes <- length(selectGenes)      
  if (is.null(dynamicParams$knots)) {
       allknots <- cbind(runif(ngenes, 0, .5), runif(ngenes, .5, 1)) # Two sites of major change
  } else {allknots <- dynamicParams$knots}
  if (is.null(dynamicParams$theta)) {
      alltheta <- matrix(rnorm(5, 5, 5), ncol=5, nrow=ngenes) # Directional Changes (num knot+degree+1)
  } else {alltheta <- dynamicParams$theta}
  
  rownames(allknots) <- selectGenes
  rownames(alltheta) <- selectGenes
  
  R1 <- sapply(selectGenes, function(x){
    genpts <- sort(runif(numCells))
    knots <- allknots[x,] 
    theta <- alltheta[x,]
    theta <- sqrt(mu[x])*scale(theta, T, F) + mu[x]
    theta[theta<0] <- 0
    basis <- bs(x = genpts, knots = knots, degree = dynamicParams$degree,
                Boundary.knots = c(0,1), intercept = TRUE)
    resp <- basis %*% theta
    use_dt <- data.table(xvals = genpts, yvals = resp[,1])
    use_dt$yvary <- rpois(rep(1,nrow(use_dt)), use_dt$yvals)
  })
  R <- t(cbind(R0,R1))
  
  cellPopulation <- rep("NotDynamic", length(mu))
  names(cellPopulation) <- names(mu)
  cellPopulation[selectGenes] <- "Dynamic"
  
  return(list(R,cellPopulation))
}