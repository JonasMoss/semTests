#' Asymptotically distribution free covariance matrix.
#' @param x Data.
#' @param sigma Covariance of the data.
#' @return Estimate of the ADF covariance matrix.
#' @keywords Internal.
gamma_est_adf <- function(x) {
  i_row <- \(n) unlist(lapply(seq_len(n), seq.int, n))
  i_col <- \(n) rep.int(seq_len(n), times = rev(seq_len(n)))
  rows <- i_row(ncol(x))
  cols <- i_col(ncol(x))
  y <- t(x) - colMeans(x, na.rm = TRUE)
  z <- y[cols, ] * y[rows, ]
  mat <- z - rowMeans(z, na.rm = TRUE)
  tcrossprod(mat) / nrow(x)
}

#' Normal theory gamma matrix
#'
#' Code obtained from `lavaan`:
#' https://github.com/yrosseel/lavaan/blob/6f047c800206d23f246d484b9522295257614222/R/lav_matrix.R
#'
#' Calculate the gamma matrix from a matrix of observations.
#' @param sigma Covariance matrix of the data.
#' @return Normal theory gamma matrix.
#' @keywords internal
gamma_est_nt <- function(sigma) {
  n <- ncol(sigma)
  lower <- lower_vec_indices(n)
  upper <- upper_vec_indices(n)
  y <- sigma %x% sigma
  out <- (y[lower, , drop = FALSE] + y[upper, , drop = FALSE]) / 2
  out[, lower, drop = FALSE] + out[, upper, drop = FALSE]
}

#' Unbiased asymptotic covariance matrix.
#'
#' @param x Data.
#' @param sigma Unbiased covariance matrix of the data.
#' @param gamma_adf The `gamma` matrix. If `NULL`, computes gamma from `x` and `sigma`.
#' @param gamma_nt The `gamma` matrix under normal theory.
#'    If `NULL`, computes `gamma_nt` from `sigma`.
#' @return Unbiased asymptotic covariance matrix.
#' @keywords internal
gamma_est_unbiased <- function(x, sigma = NULL, gamma_adf = NULL, gamma_nt = NULL) {
  n <- nrow(x)
  sigma <- (n - 1) / n * (if(is.null(sigma)) cov(x) else sigma)
  gamma_adf <- if(is.null(gamma_adf)) gamma_est_adf(x) else gamma_adf
  gamma_nt <- if(is.null(gamma_nt)) gamma_est_nt(sigma) else gamma_nt
  gamma_rem <- tcrossprod(vech(sigma))
  mult <- n / ((n - 2) * (n - 3))
  mult * ((n - 1) * gamma_adf - (gamma_nt - 2 / (n - 1) * gamma_rem))
}

#' Obtain indices of lower or upper triangular matrix in vec indices.
#'
#' Code obtained from `lavaan`:
#' https://github.com/yrosseel/lavaan/blob/6f047c800206d23f246d484b9522295257614222/R/lav_matrix.R
#'
#' @param n Dimension of square matrix.
#' @param diagonal If `TRUE`, includes the diagonal elements.
#' @returns Indices `x` so that `a[x] = c(a)[x]` returns the elements
#'    of the lower (upper) diagonal matrix in row-wise (column-wise)
#'    order.
#' @keywords internal
#'
lower_vec_indices <- function(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  if (diagonal) which(rows >= cols) else which(rows > cols)
}

upper_vec_indices <- function(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
  if (diagonal) tmp[rows >= cols] else tmp[rows > cols]
}

#' Half-vectorize matrix.
#'
#' @param x Matrix to vectorize.
#' @keywords internal
vech <- function(x) x[row(x) >= col(x)]

