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
#' @inheritParams estimateScaffoldParams
#' @import splines
generateDynamicGeneCounts <- function(numCells, mu, dynamicParams) {
  

  selectGenes <- sample(1:length(mu), ceiling(dynamicParams$propDynamic*length(mu)))
  otherGenes <- setdiff(1:length(mu), selectGenes)
  
  R0 <- matrix(sapply(otherGenes, function(x) rpois(numCells, mu[x])), nrow=numCells)

  R1 <- sapply(selectGenes, function(x){
    genpts <- sort(runif(numCells))
    knots <- c(runif(1,0, .5), runif(1,.5, 1)) # Two sites of major change
    # Directional Changes (num knot+degree+1)
    theta <- c(rnorm(1, 5, 5), rnorm(1, 5, 5), rnorm(1, 5, 5), rnorm(1, 5, 5), rnorm(1, 5, 5))
    theta <- sqrt(mu[x])*scale(theta, T, F) + mu[x]
    theta[theta<0] <- 0
    basis <- bs(x = genpts, knots = knots, degree = 2,
                Boundary.knots = c(0,1), intercept = TRUE)
    resp <- basis %*% theta
    use_dt <- data.table(xvals = genpts, yvals = resp[,1])
    use_dt$yvary <- rpois(rep(1,nrow(use_dt)), use_dt$yvals)
  })
  R <- t(cbind(R0,R1))
  return(R)
}