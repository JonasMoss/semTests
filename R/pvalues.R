#' Calculate p-values for a lavaan object.
#' @param object A `lavaan` object.
#' @export
#' @return A named vector containing the p-values `pml`. `psb`, `pfull`,
#'    `phalf`, `pcf`, `pss`.

pvalues = function(object) {

  tml <- lavaan::fitmeasures(object, "chisq")
  df <- lavaan::fitmeasures(object, "df")
  ug <- lavaan::inspect(object, "UG")
  lambdas <- Re(eigen(ug)$values[1:df])
  eigenps <- eigen_pvalues(tml, lambdas)
  data <- object@Data@X[[1]]
  colnames(data) <- object@Data@ov.names[[1]]
  pss_object <- lavaan::lavaan(object, data = data, test = "scaled.shifted")

  c(pml = unname(lavaan::fitmeasures(object, "pvalue")),
    psb = eigenps$psb,
    pfull = eigenps$pfull,
    phalf  = eigenps$phalf,
    pcf = unname(scaled_f(tml, lambdas)),
    pss = unname(lavaan::fitmeasures(pss_object, fit.measures = "pvalue.scaled")))

}

#' Calculate the scaled_f p-value.
#' @param tml Chi-square fit value from a lavaan object.
#' @param eigenvalues Eigenvalues of UG matrix.
#' @return scaled f p-value.
#' @keywords internal
scaled_f <- function(tml, eigenvalues) {
  s1 <- sum(eigenvalues)
  s2 <- sum(eigenvalues^2)
  s3 <- sum(eigenvalues^3)

  denom <- 2 * s1 * s2^2 - s1^2 * s3 + 2 * s2 * s3
  if (denom > 0) {
    d1F3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) / denom
    d2F3 <- (s1^2 * s2 + 2 * s2^2) / (s3 * s1 - s2^2) + 6
    if (d2F3 < 6) d2F3 <- Inf
    cF3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) / (s1^2 * s2 - 4 * s2^2 + 6 * s1 * s3)
  } else {
    d1F3 <- Inf
    d2F3 <- s1^2 / s2 + 4
    cF3 <- s1 * (s1^2 + 2 * s2) / (s1^2 + 4 * s2)
  }

  1 - stats::pf(tml / cF3, d1F3, d2F3)
}

#' Get eigenvalue based p-values.
#' @param tml Chi-square fit value from a lavaan object.
#' @param eigenvalues Eigenvalues of the UG matrix.
#' @return List of eigenvalue-based p-values. `psb` is Satorra--Bentler,
#'    `pfull` is based on every p-value, while `phalf` is based on the
#'    half of the (largest) p-values.
#' @keywords internal

eigen_pvalues <- function(tml, eigenvalues) {
  pfull <- CompQuadForm::imhof(tml, eigenvalues)$Qq
  psb <- CompQuadForm::imhof(tml, rep(mean(eigenvalues), length(eigenvalues)))$Qq
  k <- ceiling(length(eigenvalues) / 2)

  eigenvalues[1:k] <- mean(eigenvalues[1:k])
  eigenvalues[(k + 1):length(eigenvalues)] <- mean(eigenvalues[(k + 1):length(eigenvalues)])

  phalf <- CompQuadForm::imhof(tml, eigenvalues)$Qq
  list(pfull = pfull, phalf = phalf, psb = psb)
}
