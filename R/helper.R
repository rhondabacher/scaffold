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


