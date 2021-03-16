#' @export
redobox <- function(DATA, smallc) {

  DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
  y <- Rfast::Log(DATA)

  return(y)
}

# Regression using negative binomial. ignore zeros

#' @importFrom MASS glm.nb
quickreg.nb <- function(x,InputData)
{
  Data = InputData[[1]]
  ORIGdata = InputData[[2]]
  SeqDepth = InputData[[3]]
  X = InputData[[4]][x]

  ok <- names(which(Data[X, ] > 0)) # nonzero.

  Y <- Data[X, ]
  ok <- intersect(ok, names(Y))
  # Fit Neg. Bin. GLM
  slope <- try(glm.nb(ceiling(Y[ok]) ~ Rfast::Log(SeqDepth)[ok], maxit=10)$coef[2], silent=TRUE)
  slope <- as.numeric(slope)
  names(slope) <- X

  return(slope)
}

# Regression using negative binomial. use zeros
quickreg.nb.wz <- function(x,InputData)
{
  Data = InputData[[1]]
  ORIGdata = InputData[[2]]
  SeqDepth = InputData[[3]]
  X = InputData[[4]][x]

  ok <- names(which(Data[X, ] >= 0)) # nonzero.

  Y <- Data[X, ]
  ok <- intersect(ok, names(Y))
  # Fit Neg. Bin. GLM
  slope <- try(glm.nb(ceiling(Y[ok]) ~ Rfast::Log(SeqDepth)[ok], maxit=10)$coef[2], silent=TRUE)
  slope <- as.numeric(slope)
  names(slope) <- X

  return(slope)
}


# method to estimate parameter for capture efficiency
#' @importFrom stats dbinom optimize
#' @export
estimateCaptureEff <- function(Data, compareData, SD) {

  gdetectRate <- rowSums(Data!=0) /ncol(Data)

  Ptest <- rowMeans(Data) / nrow(Data)
  try1 <- split(Rfast::Sort(Ptest), cut(seq_along(Rfast::Sort(Ptest)), 10, labels = FALSE))
  randG <- do.call(c,lapply(1:10, function(x) sample(names(try1[[x]]), 100)))
  minFunc <- function(inGuess){
     tt <- dbinom(0, round(inGuess*nrow(Data)), Ptest[randG], log = TRUE)
     avg.detection.raw = mean(colMeans(compareData!=0))
     X = abs((1-mean(exp(tt), na.rm=T)) - avg.detection.raw)
     return(X)
   }
   simparm <- optimize(minFunc, lower=0, upper=1, tol=1e-10)$minimum
   simparm <- c(simparm, SD)
   print(simparm)

  esteff <- rnorm(ncol(Data), simparm[1], simparm[2])
  return(esteff)

}
