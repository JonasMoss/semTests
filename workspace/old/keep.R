#' Extract U from lavaan object.
#' @param object `lavaan` object.
get.U <- function(object) {
  V <- lavaan::lavTech(object, "WLS.V")[[1]]
  Delta <- lavaan:::computeDelta(object@Model)[[1]]
  E.inv <- solve(t(Delta) %*% V %*% Delta)
  V - V %*% Delta %*% E.inv %*% t(Delta) %*% V
}


## gives eigen p-value,
get.pvalue <- function(fit, usemean) { # returns p-value. can
  eigenvalues <- get.eigenvalues(fit)
  if (usemean) {
    k <- ceiling(length(eigenvalues) / 2)
    eigenvalues[1:k] <- mean(eigenvalues[1:k])
    eigenvalues[(k + 1):length(eigenvalues)] <-
      mean(eigenvalues[(k + 1):length(eigenvalues)])
  }

  CompQuadForm::imhof(fitmeasures(fit)["chisq"], eigenvalues)$Qq
}

get.pvalue.difftesting <- function(eigenvalues, Td, delta.df, usemean) {
  # strip
  eigenvalues <- eigenvalues[1:delta.df]

  if (usemean) {
    k <- ceiling(length(eigenvalues) / 2)
    eigenvalues[1:k] <- mean(eigenvalues[1:k])
    eigenvalues[(k + 1):length(eigenvalues)] <-
      mean(eigenvalues[(k + 1):length(eigenvalues)])
  }

  CompQuadForm::imhof(Td, eigenvalues)$Qq
}
