# Reduced symmetric products for positive-definite inference metrics.
# If H = R'R and Q spans ker((R Delta)'), then U = B B', B = R'Q.
# Only the trailing columns of the QR factor are formed, never the full Q.
residual_metric_factor <- function(metrics, delta, df) {
  sizes <- vapply(metrics, nrow, integer(1))
  if (sum(sizes) != nrow(delta) || df <= 0L || df > nrow(delta)) {
    stop("Residual metric and tangent dimensions do not agree.", call. = FALSE)
  }
  roots <- lapply(metrics, function(H) {
    if (nrow(H) != ncol(H) || any(!is.finite(H))) {
      stop("Inference metric must be a finite square matrix.", call. = FALSE)
    }
    # A non-PD metric does not admit this factorization. The caller retains
    # its original construction for that convention, rather than clipping H.
    tryCatch(chol((H + t(H)) / 2), error = function(e) NULL)
  })
  if (any(vapply(roots, is.null, logical(1)))) return(NULL)
  ends <- cumsum(sizes)
  indices <- Map(seq.int, ends - sizes + 1L, ends)
  A <- do.call(rbind, Map(function(R, ix) R %*% delta[ix, , drop = FALSE],
                         roots, indices))
  if (any(!is.finite(A))) stop("Model tangent is non-finite.", call. = FALSE)
  # Column scaling makes the rank decision insensitive to parameter units.
  scales <- sqrt(colSums(A^2))
  if (any(scales == 0)) stop("Model tangent is rank deficient.", call. = FALSE)
  qr_A <- qr(sweep(A, 2L, scales, "/"), tol = 1e-8)
  if (qr_A$rank != ncol(A) || nrow(A) - qr_A$rank != df) {
    stop("Model tangent rank does not agree with the test degrees of freedom.",
         call. = FALSE)
  }
  tail <- matrix(0, nrow(A), df)
  tail[cbind(qr_A$rank + seq_len(df), seq_len(df))] <- 1
  Q <- qr.qy(qr_A, tail)
  do.call(rbind, Map(function(R, ix) crossprod(R, Q[ix, , drop = FALSE]),
                     roots, indices))
}

positive_symmetric_eigenvalues <- function(x, df, context) {
  if (nrow(x) != df || ncol(x) != df || any(!is.finite(x))) {
    stop(context, " has invalid dimensions or non-finite entries.", call. = FALSE)
  }
  values <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 1e-8 * max(1, max(abs(values)))
  if (any(values < -tolerance)) {
    stop(context, " has materially negative eigenvalues.", call. = FALSE)
  }
  if (sum(values > tolerance) < df) {
    stop(context, " rank is smaller than the test degrees of freedom.", call. = FALSE)
  }
  pmax(values, 0)
}

# The reference product C^-1 S is similar to R'^-1 S R^-1 for C = R'R.
# Keep the general construction when C is not positive definite.
restriction_spectrum <- function(C, S, df, context) {
  R <- tryCatch(chol((C + t(C)) / 2), error = function(e) NULL)
  if (is.null(R)) {
    return(ugamma_eigenvalues(generalized_inverse(C) %*% S, df, context))
  }
  left <- forwardsolve(t(R), S)
  reduced <- t(forwardsolve(t(R), t(left)))
  positive_symmetric_eigenvalues(reduced, df, context)
}
