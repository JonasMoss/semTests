# Direct empirical moment projection is deliberately limited to the normal-
# likelihood complete-data ML convention. Preserve supplied NACOV, weights,
# alternative finite-sample Gamma corrections, and other first-stage methods.
can_project_raw_moments <- function(fit) {
  is_classic_nt(fit) && identical(fit@Options$likelihood, "normal") &&
    isTRUE(fit@Options$sample.cov.rescale) &&
    !isTRUE(fit@Options$gamma.n.minus.one) &&
    !isTRUE(fit@Options$gamma.unbiased) &&
    identical(fit@Data@data.type, "full") && fit@Data@nlevels == 1L &&
    !length(fit@Data@cluster) && !isTRUE(fit@Data@std.ov) &&
    !any(vapply(fit@Data@weights, length, integer(1)) > 0L) &&
    !length(fit@Data@sampling.weights) &&
    is.null(fit@call$NACOV) && is.null(fit@call$sample.cov)
}

# For centered observations x_i, the moment influence is
# z_i = (x_i, vech(x_i x_i' - S)), with the mean block optional and
# S = sum_i x_i x_i'/n. Accumulate (Z F)'(Z F)/n in bounded row tiles.
raw_moment_projection <- function(X, F, meanstructure = FALSE, tile_rows = NULL) {
  p <- ncol(X)
  n <- nrow(X)
  m <- p * (p + 1L) / 2L + if (meanstructure) p else 0L
  if (n < 2L || nrow(F) != m || !ncol(F) ||
      any(!is.finite(X)) || any(!is.finite(F))) {
    stop("Raw moment projection requires finite complete data and a conformable factor.",
         call. = FALSE)
  }
  if (is.null(tile_rows)) {
    # Budget for moment products, indexed temporaries and projected rows.
    tile_rows <- max(1L, min(512L, floor(8 * 1024^2 / (8 * (3 * m + 2 * ncol(F))))))
  }
  if (length(tile_rows) != 1L || !is.finite(tile_rows) || tile_rows < 1 ||
      tile_rows != floor(tile_rows)) {
    stop("tile_rows must be a positive integer.", call. = FALSE)
  }
  tile_rows <- min(n, tile_rows)
  X <- sweep(X, 2L, colMeans(X), "-")
  ij <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  S <- crossprod(X) / n
  center <- S[ij]
  projected_center <- if (meanstructure) c(rep(0, p), center) else center
  projected_center <- drop(crossprod(F, projected_center))
  out <- matrix(0, ncol(F), ncol(F))
  for (first in seq.int(1L, n, by = tile_rows)) {
    rows <- first:min(n, first + tile_rows - 1L)
    Xt <- X[rows, , drop = FALSE]
    Z <- Xt[, ij[, 1L], drop = FALSE] * Xt[, ij[, 2L], drop = FALSE]
    if (meanstructure) Z <- cbind(Xt, Z)
    Y <- sweep(Z %*% F, 2L, projected_center, "-")
    out <- out + crossprod(Y)
  }
  out / n
}

# Return F' Gamma_global F without constructing Gamma_global. Reuse cached
# NACOV when present; otherwise eligible ML fits need only raw-data tiles.
# The explicit source argument is internal, for parity checks and benchmarks.
projected_gamma <- function(fit, F, fallback = NULL,
                            source = c("auto", "raw", "materialized"),
                            tile_rows = NULL) {
  source <- match.arg(source)
  eligible <- can_project_raw_moments(fit)
  cached <- any(vapply(fit@SampleStats@NACOV, length, integer(1)) > 0L)
  use_raw <- source == "raw" || (source == "auto" && eligible && !cached)
  if (use_raw && !eligible) {
    stop("Direct raw moment projection is not supported for this fitting convention.",
         call. = FALSE)
  }
  weights <- as.numeric(unlist(fit@SampleStats@nobs)) / fit@SampleStats@ntotal
  blocks <- if (use_raw) {
    X <- lavaan::lavInspect(fit, "data")
    if (is.list(X)) X else list(X)
  } else {
    gamma_from_lavaan(fit, fallback)
  }
  sizes <- if (use_raw) {
    vapply(blocks, function(X) {
      p <- ncol(X)
      p * (p + 1L) / 2L + if (isTRUE(fit@Options$meanstructure)) p else 0L
    }, numeric(1))
  } else {
    vapply(blocks, nrow, integer(1))
  }
  if (length(blocks) != length(weights) || sum(sizes) != nrow(F)) {
    stop("Moment covariance blocks and projection dimensions do not agree.", call. = FALSE)
  }
  ends <- cumsum(sizes)
  indices <- Map(seq.int, ends - sizes + 1L, ends)
  out <- matrix(0, ncol(F), ncol(F))
  for (g in seq_along(blocks)) {
    Fg <- F[indices[[g]], , drop = FALSE]
    block <- if (use_raw) {
      if (nrow(blocks[[g]]) != fit@SampleStats@nobs[[g]]) {
        stop("Raw data and fitted sample sizes do not agree.", call. = FALSE)
      }
      raw_moment_projection(blocks[[g]], Fg, isTRUE(fit@Options$meanstructure), tile_rows)
    } else {
      crossprod(Fg, blocks[[g]] %*% Fg)
    }
    out <- out + block / weights[g]
  }
  out
}
