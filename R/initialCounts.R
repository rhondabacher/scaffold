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
    R <- (sapply(1:length(mu), function(x) rnegbin(numCells, mu[x], theta[x])))

    R <- t(R*Rfast::Sort(runif(nrow(R), degree[1], degree[2])))
  }
  if(type=='p') {
    R <- (sapply(1:length(mu), function(x) rpois(numCells, mu[x])))

    R <- t(R*Rfast::Sort(runif(nrow(R), degree[1], degree[2])))
  }
  return(R)
}
