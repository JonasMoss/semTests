#' Calculate p-values for a lavaan object.
#'
#' Calculate p-values for a `lavaan` object using several methods.
#'
#' * `pstd` the standard *p*-value extracted from lavaan.
#' * `psb` Satorra-Bentler *p*-value.
#' * `pfull` *p*-value based on all eigenvalues of the gamma matrix.
#' * `phalf` *p*-value based on largest half of eigenvalues of the gamma matrix.
#' * `pcf` Scaled F *p*-value.
#' * `pss` Scaled and shifted *p*-value.
#'
#' @param object A `lavaan` object.
#' @export
#' @return A named vector containing the p-values `pml`. `psb`, `pfull`,
#'    `phalf`, `pcf`, `pss`.
pvalues <- function(object) {
  chisq <- lavaan::fitmeasures(object, "chisq")
  df <- lavaan::fitmeasures(object, "df")
  ug <- lavaan::inspect(object, "UG")
  lambdas <- Re(eigen(ug)$values[1:df])
  eigenps <- eigen_pvalues(chisq, lambdas)
  c(
    pstd = unname(lavaan::fitmeasures(object, "pvalue")),
    psb = eigenps$psb,
    pfull = eigenps$pfull,
    phalf = eigenps$phalf,
    plog = eigenps$plog,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(object),
    pmv = mean_var_adjusted(object)
  )
}

#' Calculate the scaled and shifted / the mean-variance adjusted p-value
#'
#' Code copied from `lavaan:::lav_test_satorra_bentler`.
#'
#' @param object `lavaan` object.
#' @name laavan_tests
#' @return The scaled and shifted p-value or the mean-variance adjusted p-value.
NULL

#' @rdname laavan_tests
scaled_and_shifted <- function(object) {

  group <- object@test$standard$stat.group
  lavsamplestats <- object@SampleStats
  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
  df <- object@test$standard$df
  ug <- lavaan::inspect(object, "UG")

  trace_ug <- sum(diag(ug))
  trace_ug2 <- sum(diag(ug %*% ug))
  a <- sqrt(df / trace_ug2)
  shift_parameter <- fg * (df - a * trace_ug)
  scaling_factor <- 1 / a
  if (scaling_factor < 0) scaling_factor <- as.numeric(NA)
  stat_group <- (group * a + shift_parameter)
  stat <- sum(stat_group)
  unname(1 - stats::pchisq(stat, df))

}

#' @rdname laavan_tests
mean_var_adjusted <- function(object) {
  ug <- lavaan::inspect(object, "UG")
  group <- object@test$standard$stat.group

  trace_ug <- sum(diag(ug))
  trace_ug2 <- sum(diag(ug %*% ug))
  df <- trace_ug^2 / trace_ug2
  scaling_factor <- trace_ug / df
  if (scaling_factor < 0) scaling_factor <- as.numeric(NA)
  stat_group <- group / scaling_factor
  stat <- sum(stat_group)
  unname(1 - stats::pchisq(stat, df))
}

#' Calculate the scaled_f p-value.
#' @param chisq Chi-square fit value from a lavaan object.
#' @param eig eig of UG matrix.
#' @return scaled f p-value.
#' @keywords internal
scaled_f <- function(chisq, eig) {
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

#' Get eigenvalue-based p-values.
#' @param chisq Chi-square fit value from a lavaan object.
#' @param eig Eigenvalues of the UG matrix.
#' @return List of eigenvalue-based p-values. `psb` is Satorra--Bentler,
#'    `pfull` is based on every p-value, while `phalf` is based on the
#'    half of the (largest) p-values.
#' @keywords internal
eigen_pvalues <- function(chisq, eig) {
  m <- length(eig)
  pfull <- CompQuadForm::imhof(chisq, eig)$Qq
  psb <- CompQuadForm::imhof(chisq, rep(mean(eig), m))$Qq
  alpha <- 0.0
  lambdas_log <- mean(eig) + (1 - alpha) * (mean(log(seq(m))) - log(seq(m)))
  plog <- CompQuadForm::imhof(chisq, lambdas_log)$Qq

  k <- ceiling(m / 2)
  eig[1:k] <- mean(eig[1:k])
  eig[(k + 1):m] <- mean(eig[(k + 1):m])

  phalf <- CompQuadForm::imhof(chisq, eig)$Qq
  list(pfull = pfull, phalf = phalf, psb = psb, plog = plog)
}
