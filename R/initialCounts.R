#' Generate initial counts
#'
#' Generates the transcripts in the cell. 
#'
#' @param numCells The number of cells to generate counts for.
#' @param mu a vector containing the mean of each gene to simulate counts for.
#' @param theta only needed if using the negative binomial (type='nb'). This is the theta
#' where the variance is mu + mu^2/theta
#' @param type the underlying distribution to simulate counts from. 
#' Highly recomend sticking with the default Poisson, type='p'.
#' @param degree the degree of heterogeneity within the cell population.
#'	Should be a vector of length two specifying the bounds
#'
#' @return Output is matrix of counts, genes on the rows and cells
#'   on the columns.
#'
#' @importFrom stats rpois
#' @importFrom MASS rnegbin
#' @export
#' @examples 
#'  # generate 50 cells each having 1000 genes
#' mu <- rnorm(1000, 100, 25)
#' initmat <- generateGeneCounts(50, mu = mu, type='p', degree=c(1/2,2))
#' head(initmat)
#'
#' @export

generateGeneCounts <- function(numCells, mu, theta, type, degree) {
	
 # multiple options just in case, mainly use the poisson option.
  if(type=='nb') {
    R <- matrix(sapply(1:length(mu), function(x) rnegbin(numCells, mu[x], theta[x])), nrow=numCells)

    R <- t(R*Rfast::Sort(runif(nrow(R), degree[1], degree[2])))
  }
  if(type=='p') {
    R <- matrix(sapply(1:length(mu), function(x) rpois(numCells, mu[x])), nrow=numCells)

    R <- t(R*Rfast::Sort(runif(nrow(R), degree[1], degree[2])))
  }
  return(R)
}

## The motivation for this code is from the simstudy package.
generateDynamicGeneCounts <- function(numCells, mu, propDynamic) {
  
	library(splines)
  selectGenes <- sample(1:length(mu), ceiling(propDynamic*length(mu)))
  otherGenes <- setdiff(1:length(mu), selectGenes)
  
  R0 <- matrix(sapply(otherGenes, function(x) rpois(numCells, mu[x])), nrow=numCells)
  dim(R0)
  length(otherGenes)
  
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