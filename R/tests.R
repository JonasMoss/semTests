#' Calculate the scaled and shifted / the mean-variance adjusted p-value
#'
#' @param chisq Chi-square fit value from a lavaan object.
#' @param lambdas Eigenvalues of UG matrix.
#' @name laavan_tests
#' @keywords internal
#' @return The scaled and shifted p-value or the mean-variance adjusted p-value.
NULL

#' @rdname laavan_tests
scaled_and_shifted <- \(chisq, lambdas) {
  df <- length(lambdas)
  tr_ug <- sum(lambdas)
  tr_ug2 <- sum(lambdas^2)
  a <- sqrt(df / tr_ug2)
  b <- df - sqrt(df * tr_ug^2 / tr_ug2)
  t3 <- unname(chisq * a + b)
  1 - stats::pchisq(t3, df = df)
}

#' Calculate the scaled_f p-value.
#' @param chisq Chi-square fit value from a lavaan object.
#' @param eig eig of UG matrix.
#' @return scaled f p-value.
#' @keywords internal
scaled_f <- \(chisq, eig) {
  s1 <- sum(eig)
  s2 <- sum(eig^2)
  s3 <- sum(eig^3)
  denom <- 2 * s1 * s2^2 - s1^2 * s3 + 2 * s2 * s3
  if (denom > 0) {
    d1f3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) / denom
    d2f3 <- (s1^2 * s2 + 2 * s2^2) / (s3 * s1 - s2^2) + 6
    if (d2f3 < 6) d2f3 <- Inf
    cf3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) /
      (s1^2 * s2 - 4 * s2^2 + 6 * s1 * s3)
  } else {
    d1f3 <- Inf
    d2f3 <- s1^2 / s2 + 4
    cf3 <- s1 * (s1^2 + 2 * s2) / (s1^2 + 4 * s2)
  }
  unname(1 - stats::pf(chisq / cf3, d1f3, d2f3))
}

#' Calculate the jth eba pvalue.
#' @keywords internal
eba_pvalue <- \(chisq, lambdas, j) {
  m <- length(lambdas)
  k <- ceiling(m / j)
  eig <- lambdas
  eig <- c(eig, rep(NA, k * j - length(eig)))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  repeated <- rep(eig_means, each = k)[seq(m)]
  CompQuadForm::imhof(chisq, repeated)$Qq
}

#' Calculate the jth all pvalue.
#' @keywords internal
pvalue_all <- \(chisq, lambdas) {
  CompQuadForm::imhof(chisq, lambdas)$Qq
}

#' Calculate the jth pall pvalue.
#' @keywords internal
pall <- \(chisq, lambdas) {
  CompQuadForm::imhof(chisq, lambdas/2 + mean(lambdas)/2)$Qq
}

#' Calculate the jth eba pvalue.
#' @keywords internal
peba_pvalue <- \(chisq, lambdas, j) {
  m <- length(lambdas)
  k <- ceiling(m / j)
  eig <- lambdas
  eig <- c(eig, rep(NA, k * j - length(eig)))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  eig_mean <- mean(lambdas)
  repeated <- rep(eig_means, each = k)[seq(m)]
  CompQuadForm::imhof(chisq, (repeated + eig_mean) / 2)$Qq
}

#' Calculate penalized OLS pvalue.
#' @keywords internal
pols_pvalue <- \(chisq, lambdas, gamma) {

  lambda_hat <- if(length(lambdas) == 1) {
    lambdas
  } else {
    x <- seq_along(lambdas)
    beta1_hat <- 1 / gamma * stats::cov(x, lambdas) / stats::var(x)
    beta0_hat <- mean(lambdas) - beta1_hat * mean(x)
    lambda_hat <- pmax(beta0_hat + beta1_hat * x, 0)
  }
  CompQuadForm::imhof(chisq, lambda_hat)$Qq
}
