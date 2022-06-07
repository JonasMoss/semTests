#' Turn eigenvalus into blocks.
#'
#' @param eigenvalues Vector of eigenvalues.
#' @param blocks Number of desired blocks. Must be greater than 1.
#' @return Blocked eigenvalues.

eigen.block <- function(eigenvalues, blocks) {

  stopifnot(blocks > 1)

  e <- eigenvalues
  t <- chunk(eigenvalues, blocks)
  m <- sapply(t, mean)

  counter <- 1

  for (i in 1:length(t)) {
    len <- length(t[[i]])
    e[counter:(counter + len - 1)] <- m[i]
    counter <- counter + len
  }

  e

}


#' Difference of eigenvalue test
#'
#' @param object_0,object_1 `lavaan` objects to compare.
#' @param fit `lavaan` object.
#' @param use_mean Use the mean?
#' @return Unsure.

get.eigen.diff <- function(object_0, object_1, use_mean) { # returns eigenvalues of U_d Gamma
  ud_gamma <- inspect(object_0, "UGamma") - inspect(object_1, "UGamma")

  df <- fitMeasures(object_0, "df") - fitMeasures(object_1, "df")
  eigenvalues <- Re(eigen(ud_gamma)$values[1:df]) # only the first. BUT: there remaining are not zero?

  if (use_mean) {
    k <- ceiling(length(eigenvalues) / 2)
    eigenvalues[1:k] <- mean(eigenvalues[1:k])
    eigenvalues[(k + 1):length(eigenvalues)] <-
      mean(eigenvalues[(k + 1):length(eigenvalues)])
  }

  CompQuadForm::imhof(lavaan::fitmeasures(fit)["chisq"], eigenvalues)$Qq

}


eigen.jenks <- function(eigenvalues, numblocks) {

  e <- eigenvalues
  max.e <- max(e)
  t <- BAMMtools::getJenksBreaks(eigenvalues, k = numblocks + 1)

  for (i in 2:length(t)) {
    indx <- which(eigenvalues <= t[i])
    e[indx] <- mean(eigenvalues[indx])
    eigenvalues[indx] <- max.e + 1
  }

  e

}

eigenvalues_diff <- function(m1, m0, A.method = "exact") {
  # or delta. Note that shell command lavTestLRT has exact, while lav_test_diff_Satorra2000 has default delta.

  # extract information from m1 and m2
  T1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df

  T0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df

  # m = difference between the df's'
  m <- r0 - r1

  Gamma <- lavTech(m1, "Gamma") # the same for m1 and m0

  WLS.V <- lavTech(m1, "WLS.V")
  PI <- lavaan:::computeDelta(m1@Model)
  P <- lavTech(m1, "information")
  # needed? (yes, if H1 already has eq constraints)
  P.inv <- lavaan:::lav_model_information_augment_invert(m1@Model,
    information = P,
    inverted = TRUE
  )
  if (inherits(P.inv, "try-error")) {
    cat("Error! in P.inv \n")
    return(NA)
  }

  A <- lavaan:::lav_test_diff_A(m1, m0, method = A.method, reference = "H1")


  APA <- A %*% P.inv %*% t(A)
  cSums <- colSums(APA)
  rSums <- rowSums(APA)
  empty.idx <- which(abs(cSums) < .Machine$double.eps^0.5 &
    abs(rSums) < .Machine$double.eps^0.5)
  if (length(empty.idx) > 0) {
    A <- A[-empty.idx, , drop = FALSE]
  }

  # PAAPAAP
  PAAPAAP <- P.inv %*% t(A) %*% solve(A %*% P.inv %*% t(A)) %*% A %*% P.inv

  g <- 1
  UG.group <- WLS.V[[g]] %*% Gamma[[g]] %*% WLS.V[[g]] %*%
    PI[[g]] %*% PAAPAAP %*% t(PI[[g]])


  return(Re(eigen(UG.group)$values))
}

eigen.dwls = function(my.sample, my.model){
  pstar = ncol(my.sample)*(ncol(my.sample)+1)/2L
  W=matrix(rep(0, pstar*pstar), ncol=pstar)
  diag(W) =  1/diag(lavaan:::lav_samplestats_Gamma(my.sample))
  f=sem(bollen.model, data=my.sample, estimator="WLS", WLS.V=W)
  UG = inspect(f, "UGamma")
  return(Re(eigen(UG)$values[1:35]))
}
