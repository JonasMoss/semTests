#' Generate list of Vale--Maurelli samples.
#' @param n Number of samples, defaults to `n = 100`.
#' @param reps Number of repetitions.
#' @param cov Covariance matrix.
#' @param skewness Marginal skewness.
#' @param kurtosis Marginal kurtosis.
#' @return List of `n` Vale--Maurelli samples.

vale_maurelli <- function(n = 100, reps, cov, skewness, kurtosis) {

  cor = cov2cor(cov)

  list.of.samples <- list() # to be returned

  if (skewness==0 & kurtosis==0){

    for (i in 1:reps){
      X <- data.frame(mvrnorm(n, mu=rep(0,ncol(cor)), Sigma=cor))
      colnames(X) <- colnames(cor)
      list.of.samples[[length(list.of.samples)+1]] <-
        scale(X, center = FALSE, scale = 1/sqrt(diag(cov)))[,,drop=FALSE]
    }

    return(list.of.samples)

  }

  #otherwise, not normal:
  fleishman1978_abcd <- function(skewness, kurtosis) {

    system.function <- function(x, skewness, kurtosis) {
      b.=x[1L]; c.=x[2L]; d.=x[3L]
      eq1 <- b.*b. + 6*b.*d. + 2*c.*c. + 15*d.*d. - 1
      eq2 <- 2*c.*(b.*b. + 24*b.*d. + 105*d.*d. + 2) - skewness
      eq3 <- 24*(b.*d. + c.*c.*(1 + b.*b. + 28*b.*d.) +
                   d.*d.*(12 + 48*b.*d. + 141*c.*c. + 225*d.*d.)) - kurtosis
      eq <- c(eq1,eq2,eq3)
      sum(eq*eq) ## SS
    }

    out <- nlminb(start=c(1,0,0), objective=system.function,
                  scale=10,
                  control=list(trace=0),
                  skewness=skewness, kurtosis=kurtosis)
    if(out$convergence != 0) warning("no convergence")
    b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
    c(a.,b.,c.,d.)

  }

  getIcov <- function(b1, c1, d1, b2, c2, d2, R) {

    objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {
      rho=x[1L]
      eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) +
        rho*rho*(2*c1*c2) + rho*rho*rho*(6*d1*d2) - R
      eq*eq
    }

    out <- nlminb(start=R, objective=objectiveFunction,
                  scale=10, control=list(trace=0),
                  b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)
    if(out$convergence != 0) warning("no convergence")
    rho <- out$par[1L]
    rho

  }

  # number of variables
  nvar <- ncol(cor)
  # check skewness
  if(is.null(skewness)) {
    SK <- rep(0, nvar)
  } else if(length(skewness) == nvar) {
    SK <- skewness
  } else if(length(skewness) == 1L) {
    SK <- rep(skewness, nvar)
  } else {
    stop("skewness has wrong length")
  }

  if(is.null(kurtosis)) {
    KU <- rep(0, nvar)
  } else if(length(kurtosis) == nvar) {
    KU <- kurtosis
  } else if(length(kurtosis) == 1L) {
    KU <- rep(kurtosis, nvar)
  } else {
    stop("kurtosis has wrong length")
  }

  # create Fleishman table
  FTable <- matrix(0, nvar, 4L)
  for(i in 1:nvar) {
    FTable[i,] <- fleishman1978_abcd(skewness=SK[i], kurtosis=KU[i])
  }

  # compute intermediate correlations between all pairs
  Icor <- diag(nvar)
  for(j in 1:(nvar-1L)) {
    for(i in (j+1):nvar) {
      if(cor[i,j] == 0) next
      Icor[i,j] <- Icor[j,i] <-
        getIcov(FTable[i,2], FTable[i,3], FTable[i,4],
                FTable[j,2], FTable[j,3], FTable[j,4], R=cor[i,j])
    }
  }

  list.of.samples <- list()
  for (i in 1:reps){

    # generate Z ## FIXME: replace by rmvnorm once we use that package
    X <- Z <- MASS::mvrnorm(n=n, mu=rep(0,nvar), Sigma=Icor)

    # transform Z using Fleishman constants
    for(i in 1:nvar) {
      X[,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]*Z[,i] +
        FTable[i,4L]*Z[,i]*Z[,i]*Z[,i]
    }

    X= data.frame(X)
    colnames(X) <- colnames(cor)

    #scale
    # first, we scale


    list.of.samples[[length(list.of.samples)+1]] <- scale(X, center = FALSE,
                                                          scale = 1/sqrt(diag(cov)))[,,drop=FALSE]

  }

  return(list.of.samples)
}

