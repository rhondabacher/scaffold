#' @importFrom stats rpois
#' @importFrom MASS rnegbin
generateGeneCounts <- function(numCells, mu, theta, type, degree) {

 #multiple options just in case, mainly use the poisson option.
  if(type=='nb') {
    R <- (sapply(1:length(mu), function(x) rnegbin(numCells, mu[x], theta[x])))

    R <- t(R*sort(runif(nrow(R), 1/degree, degree)))
  }
  if(type=='p') {
    R <- (sapply(1:length(mu), function(x) rpois(numCells, mu[x])))

    R <- t(R*sort(runif(nrow(R), 1/degree, degree)))
  }
  if(type=='np') {
    R <- c()
    sj <- sort(runif(numCells, 1/degree, degree))

    for(j in 1:length(sj)) {
      R <- cbind(R, sapply(1:length(mu), function(x) sj[j]*rnegbin(1, mu[x], theta[x])))
    }
  }

  return(R)
}
