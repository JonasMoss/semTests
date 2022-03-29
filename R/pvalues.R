#' Refitted p-values for a lavaan object.
#'
#' @param object A fitted `lavaan` object.
#' @param pvalue Passed to `lavaan::fitmeasures` as `fit.measures`.
#'   Defaults to `pvalue.scaled`.
#' @param ... Passed to `lavaan::sem`.
refitted_pvalue <- function(
  object,
  pvalue = c("pvalue.scaled", "pvalue"),
  ...) {
  pvalue <- match.arg(pvalue)
  data <- object@Data@X[[1]]
  colnames(data) <- object@Data@ov.names[[1]]
  object_ <- lavaan::sem(object, data = data, ...)
  unname(lavaan::fitmeasures(object_, fit.measures = pvalue))
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
scaled_and_shifted = function(object) {
  df <- object@test$standard$df
  ug <- lavaan::inspect(object, "UG")
  trace_ug <- sum(diag(ug))
  trace_ug2 <- sum(diag(ug %*% ug))
  group <- object@test$standard$stat.group
  lavsamplestats = object@SampleStats
  fg <- unlist(lavsamplestats@nobs)/lavsamplestats@ntotal
  a <- sqrt(df/trace_ug2)
  shift_parameter <- fg * (df - a*trace_ug)
  scaling_factor  <- 1/a
  if(scaling_factor < 0) scaling_factor <- as.numeric(NA)
  stat.group <- (group * a + shift_parameter)
  stat <- sum(stat.group)
  unname(1 - stats::pchisq(stat, df))
}

#' @rdname laavan_tests
mean_var_adjusted = function(object) {
  ug <- lavaan::inspect(object, "UG")
  trace_ug <- sum(diag(ug))
  trace_ug2 <- sum(diag(ug %*% ug))
  df <- trace_ug ^ 2 / trace_ug2
  scaling_factor <- trace_ug/df
  group <- object@test$standard$stat.group
  if(scaling_factor < 0) scaling_factor <- as.numeric(NA)
  stat_group <- group/ scaling_factor
  stat <- sum(stat_group)
  unname(1 - stats::pchisq(stat, df))
}

#' Calculate p-values for a lavaan object.
#'
#' Calculate p-values for a `lavaan` object using several methods.
#'
#' * `pml` the standard *p*-value extracted from lavaan.
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
  tml <- lavaan::fitmeasures(object, "chisq")
  df <- lavaan::fitmeasures(object, "df")
  ug <- lavaan::inspect(object, "UG")
  lambdas <- Re(eigen(ug)$values[1:df])
  eigenps <- eigen_pvalues(tml, lambdas)
  c(
    pml = unname(lavaan::fitmeasures(object, "pvalue")),
    psb = eigenps$psb,
    pfull = eigenps$pfull,
    phalf = eigenps$phalf,
    plog = eigenps$plog,
    pcf = scaled_f(tml, lambdas),
    pss = scaled_and_shifted(object),
    pmv = mean_var_adjusted (object)
  )
}

#' Calculate the scaled_f p-value.
#' @param tml Chi-square fit value from a lavaan object.
#' @param eig eig of UG matrix.
#' @return scaled f p-value.
#' @keywords internal
scaled_f <- function(tml, eig) {
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

  unname(1 - stats::pf(tml / cf3, d1f3, d2f3))
}

#' Get eigenvalue-based p-values.
#' @param tml Chi-square fit value from a lavaan object.
#' @param eig Eigenvalues of the UG matrix.
#' @return List of eigenvalue-based p-values. `psb` is Satorra--Bentler,
#'    `pfull` is based on every p-value, while `phalf` is based on the
#'    half of the (largest) p-values.
#' @keywords internal
eigen_pvalues <- function(tml, eig) {
  m <- length(eig)
  pfull <- CompQuadForm::imhof(tml, eig)$Qq
  psb <- CompQuadForm::imhof(tml, rep(mean(eig), m))$Qq
  alpha <- 0.0
  lambdas_log <- mean(eig) + (1 - alpha)* (mean(log(seq(m))) - log(seq(m)))
  plog <- CompQuadForm::imhof(tml, lambdas_log)$Qq

  k <- ceiling(m / 2)
  eig[1:k] <- mean(eig[1:k])
  eig[(k + 1):m] <- mean(eig[(k + 1):m])

  phalf <- CompQuadForm::imhof(tml, eig)$Qq
  list(pfull = pfull, phalf = phalf, psb = psb, plog = plog)
}
