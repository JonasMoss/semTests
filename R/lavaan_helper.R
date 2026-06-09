# Numerical helpers derived from lavaan (https://github.com/yrosseel/lavaan).
# The eigenvalue and FIML machinery needs a Moore-Penrose generalized inverse
# and an orthogonal-complement basis; both are used by gamma.R, get_a_matrix.R,
# and fiml_fmg.R.


generalized_inverse <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x))) {
    stop("'x' must be a numeric or complex matrix")
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  x_svd <- svd(x)
  if (is.complex(x)) {
    x_svd$u <- Conj(x_svd$u)
  }
  positive <- x_svd$d > max(tol * x_svd$d[1L], 0)
  if (all(positive)) {
    x_svd$v %*% (1 / x_svd$d * t(x_svd$u))
  } else if (!any(positive)) {
    array(0, dim(x)[2L:1L])
  } else {
    x_svd$v[, positive, drop = FALSE] %*%
      ((1 / x_svd$d[positive]) * t(x_svd$u[, positive, drop = FALSE]))
  }
}

get_orthogonal_complement <- function(mat) {
  qr_decomp <- qr(mat)
  q_mat <- qr.Q(qr_decomp, complete = TRUE)
  q_mat[, -seq_len(qr_decomp$rank), drop = FALSE]
}
