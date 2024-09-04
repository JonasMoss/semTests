#' Asymptotically distribution free covariance matrix.
#' @param x Data.
#' @return Estimate of the ADF covariance matrix.
#' @keywords internal
gamma_est_adf <- \(x) {
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
gamma_est_nt <- \(sigma) {
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
gamma_est_unbiased <- \(x, n = NULL, sigma = NULL, gamma_adf = NULL, gamma_nt = NULL) {
  if (!missing(x)) n <- nrow(x)
  sigma <- if (is.null(sigma)) stats::cov(x) * (n - 1) / n else sigma
  gamma_adf <- if (is.null(gamma_adf)) gamma_est_adf(x) else gamma_adf
  gamma_nt <- if (is.null(gamma_nt)) gamma_est_nt(sigma) else gamma_nt
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
lower_vec_indices <- \(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  if (diagonal) which(rows >= cols) else which(rows > cols)
}

upper_vec_indices <- \(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
  if (diagonal) tmp[rows >= cols] else tmp[rows > cols]
}

#' Half-vectorize matrix.
#'
#' @param x Matrix to vectorize.
#' @keywords internal
vech <- \(x) x[row(x) >= col(x)]

#' Calculate unbiased gamma from lavaan object and a gamma matrix.
#'
#' WORKS ONLY FOR MODELS WITH NO MEAN STRUCTURE.
#' @param obj,gamma Object and gamma matrix.
#' @return Unbiased gamma.
#' @keywords internal
gamma_unbiased <- \(obj, gamma) {
  gamma_est_unbiased(
    n = lavaan::lavInspect(obj, "nobs"),
    sigma = obj@SampleStats@cov[[1]],
    gamma_adf = gamma,
    gamma_nt = NULL
  )
}

#' Calculate nested ugamma.
#'
#' This can also be used with restrictions.
#'
#' @param m0,m1 Two nested `lavaan` objects.
#' @param a The `A` matrix. If if `NULL`, gets calculated by
#'    `lavaan:::lav_test_diff_A` with `method = method`.
#' @param method Method passed to `lavaan:::lav_test_diff_A`.
#' @param unbiased If `TRUE`, uses the unbiased gamma estimate.
#' @keywords internal
#' @return Ugamma for nested object.
lav_ugamma_nested_2000 <- \(m0, m1, gamma, a = NULL, method = "delta", unbiased = FALSE) {
  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")

  p_inv <- lavaan::lavInspect(m1, what = "inverted.information")

  if (is.null(a)) {
    a <- lavaan:::lav_test_diff_A(m1, m0, method = method, reference = "H1")
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  paapaap <- p_inv %*% t(a) %*% MASS::ginv(a %*% p_inv %*% t(a)) %*% a %*% p_inv

  # Compute scaling factor
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

  # We need the global gamma, cf. eq.~(10) in Satorra (2000).
  # And the global V.
  # Note that the calculation is not the same as in lavaan(2000), as we
  # need more than the trace to be correct. (The simplications used in lavaan
  # do not hold for the whole ugamma.)

  for (i in (seq_along(gamma))) {
    gamma[[i]] <- 1/fg[i] * gamma[[i]]
    wls_v[[i]] <- fg[i] * wls_v[[i]]
  }

  gamma_global <- as.matrix(Matrix::bdiag(gamma))
  v_global <- as.matrix(Matrix::bdiag(wls_v))
  pi_global <- if (is.list(pi)) {
    do.call(rbind, pi)
  } else {
    pi
  }

  # U global version, eq.~(22) in Satorra (2000).
  u_global <- v_global %*% pi_global %*% paapaap %*% t(pi_global) %*% v_global
  return(u_global %*% gamma_global)
}

get_gamma <- \(m1, unbiased = FALSE, collapse = TRUE, m0 = NULL) {
  stopifnot(is.logical(unbiased))

  get_subgamma <- \(object) {
    lavoptions = lavaan::lavInspect(object, "options")
    lavdata = object@Data
    lavoptions$gamma.unbiased <- unbiased
    lavaan:::lav_samplestats_from_data(lavdata, lavoptions)@NACOV
  }

  gamma_list <- get_subgamma(m1)
  if(is.null(unlist(gamma_list))) {
    if (!is.null(m0)) {
      gamma_list <- get_subgamma(m0)
      if (is.null(unlist(gamma_list))) {
        stop("Could not calculate the gamma matrix from `m0` or `m1`. Use either `estimator = \"MLM\"' or `test=\"satorra.bentler\"' when fitting your lavaan model.")
      }
    } else {
        stop("Could not calculate the gamma matrix from the lavaan object. Use either `estimator = \"MLM\"' or `test=\"satorra.bentler\"' when fitting your lavaan model.")
    }
  }


  if (collapse) {
    fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal
    gamma_rescaled <- gamma_list
    for (i in (seq_along(gamma_list))) {
      gamma_rescaled[[i]] <- 1 / fg[i] * gamma_rescaled[[i]]
    }
    lavaan:::lav_matrix_bdiag(gamma_rescaled)
  } else {
    gamma_list
  }
}

gamma_unbiased <- \(obj, gamma) {
  gamma_est_unbiased(
    n = lavaan::lavInspect(obj, "nobs"),
    sigma = obj@SampleStats@cov[[1]],
    gamma_adf = gamma,
    gamma_nt = NULL
  )
}

