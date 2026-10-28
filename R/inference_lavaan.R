# Complete-data normal-likelihood adapters. Keep estimator and hypothesis
# boundaries explicit while missing-data and misspecification adapters develop.
check_ml_inference_adapter <- function(object, context, allow_fiml = FALSE) {
  check_supported(object)
  raw_ok <- can_project_raw_moments(object)
  if (allow_fiml && is_fiml(object)) {
    raw_ok <- identical(object@Options$estimator, "ML") &&
      identical(object@Options$likelihood, "normal") &&
      identical(object@Data@data.type, "full") && object@Data@nlevels == 1L &&
      !length(object@Data@cluster) && !isTRUE(object@Data@std.ov) &&
      !any(lengths(object@Data@weights) > 0L) && !length(object@Data@sampling.weights) &&
      is.null(object@call$sample.cov) && is.null(object@call$NACOV)
  }
  if (isTRUE(object@Options$group.w.free)) raw_ok <- FALSE
  if (!raw_ok) {
    semtests_abort(paste0(context, " currently requires complete-data normal-likelihood ML ",
      "with raw data, joint random-x inference, standard covariance normalization, ",
      "and no sampling weights, clusters or supplied sample covariance/NACOV."),
      "semTests_error_unsupported_inference")
  }
  if (length(object@Model@ceq.nonlinear.idx) || nrow(object@Model@cin.JAC) > 0L ||
      any(object@ParTable$op %in% c("<", ">")) ||
      any(is.finite(object@ParTable$lower)) || any(is.finite(object@ParTable$upper))) {
    semtests_abort(paste0(context, " currently supports affine equality constraints only; ",
                         "inequalities, bounds and nonlinear constraints are not supported."),
                  "semTests_error_unsupported_inference")
  }
}

# Fisher information in (mean, vech(covariance)) coordinates. Covariance-score
# weights are 1/2 on the diagonal and 1 off diagonal; no Kronecker matrix.
gaussian_moment_geometry <- function(Sigma, meanstructure) {
  p <- nrow(Sigma)
  K <- chol2inv(inference_chol(as.matrix(Sigma), "Implied covariance"))
  ij <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  w <- ifelse(ij[, 1L] == ij[, 2L], .5, 1)
  i <- ij[, 1L]
  j <- ij[, 2L]
  Hcov <- (K[i, i, drop = FALSE] * K[j, j, drop = FALSE] +
           K[i, j, drop = FALSE] * K[j, i, drop = FALSE]) * tcrossprod(w)
  H <- if (meanstructure) lavaan::lav_matrix_bdiag(K, Hcov) else Hcov
  list(precision = K, index = ij, weight = w, information = H,
       meanstructure = meanstructure)
}

# X is centered at the fitted mean, or at the sample mean for covariance-only
# fitting. Centering at the sample mean when the fitted mean differs is wrong.
gaussian_score_rows <- function(X, geometry) {
  q <- X %*% geometry$precision
  ij <- geometry$index
  cov_scores <- q[, ij[, 1L], drop = FALSE] * q[, ij[, 2L], drop = FALSE]
  cov_scores <- sweep(cov_scores, 2L, geometry$precision[ij], "-")
  cov_scores <- sweep(cov_scores, 2L, geometry$weight, "*")
  if (geometry$meanstructure) cbind(q, cov_scores) else cov_scores
}

ml_global_score_parts <- function(object, score_covariance, tile_rows = 256L) {
  X <- lavaan::lavInspect(object, "data")
  if (!is.list(X)) X <- list(X)
  moments <- lavaan::lavInspect(object, "implied")
  if (!is.null(moments$cov)) moments <- list(moments)
  means <- isTRUE(object@Options$meanstructure)
  n <- as.numeric(unlist(object@SampleStats@nobs))
  N <- sum(n)
  delta <- model_effective_delta(object)
  df <- as.integer(fit_df(object))
  if (df < 1L) {
    semtests_abort("The fitted model has no testable saturated complement.",
                  "semTests_error_invalid_inference")
  }
  geometry <- lapply(moments, function(m) gaussian_moment_geometry(m$cov, means))
  H <- Map(function(g, count) g$information * count, geometry, n)
  # B B' is the residual metric; G = H^-1 B supplies score coordinates with
  # G'HG = I and Delta'HG = 0, including cross-group equality constraints.
  B <- residual_metric_factor(H, delta, df)
  if (is.null(B)) {
    semtests_abort("Saturated Fisher information must be positive definite.",
                  "semTests_error_unstable_spectrum")
  }
  sizes <- vapply(H, nrow, integer(1))
  ends <- cumsum(sizes)
  indices <- Map(seq.int, ends - sizes + 1L, ends)
  G <- do.call(rbind, Map(function(Hg, ix) {
    root <- chol(Hg)
    backsolve(root, forwardsolve(t(root), B[ix, , drop = FALSE]))
  }, H, indices))
  score <- numeric(df)
  score_cov <- matrix(0, df, df)
  full_score <- numeric(nrow(delta))
  for (g in seq_along(X)) {
    if (nrow(X[[g]]) != n[g] || any(!is.finite(X[[g]])) ||
        !identical(colnames(X[[g]]), colnames(moments[[g]]$cov))) {
      semtests_abort("Raw observations and implied-moment ordering do not agree.",
                    "semTests_error_invalid_inference")
    }
    center <- if (means) as.numeric(moments[[g]]$mean) else colMeans(X[[g]])
    Xg <- sweep(X[[g]], 2L, center, "-")
    ix <- indices[[g]]
    Gg <- G[ix, , drop = FALSE]
    group_score <- numeric(df)
    group_cov <- matrix(0, df, df)
    for (first in seq.int(1L, n[g], by = tile_rows)) {
      rows <- first:min(n[g], first + tile_rows - 1L)
      S <- gaussian_score_rows(Xg[rows, , drop = FALSE], geometry[[g]])
      Y <- S %*% Gg
      full_score[ix] <- full_score[ix] + colSums(S)
      group_score <- group_score + colSums(Y)
      group_cov <- group_cov + crossprod(Y)
    }
    if (score_covariance == "centered") {
      group_cov <- group_cov - tcrossprod(group_score) / n[g]
    }
    score <- score + group_score
    score_cov <- score_cov + group_cov
  }
  list(score = score, metric = diag(df), score_covariance = (score_cov + t(score_cov)) / 2,
       nobs = N, projection = G, tangent = delta, information = H,
       nuisance_gradient = drop(crossprod(delta, full_score)) / N,
       group_nobs = n)
}

#' Global score test for an ML or FIML lavaan model
#'
#' Test exact fit against saturated means and covariances. Gaussian likelihood
#' scores include the fitted-mean residual and are projected over the model
#' tangent. FIML contributions use each observed missingness pattern, with
#' scores projected before their cross-products are accumulated.
#'
#' By default the quadratic uses expected saturated Fisher information. `curvature`
#' selects expected or observed *saturated* curvature for nuisance projection,
#' matching magmaan's global-score convention. This is an exact-fit null.
#' For full model curvature and restricted pseudo-target inference under fixed
#' misspecification, use [score_test_nested()].
#'
#' Supports normal-likelihood ML/MLM/MLR and FIML, one or several groups,
#' affine equalities, unweighted raw data, and joint random-x inference.
#' Categorical data, sampling weights, clusters, bounds and nonlinear
#' constraints are excluded. A fully studentized score uses a chi-square
#' reference. Other metrics use the full eigenvalue mixture by default.
#'
#' @param object A converged `lavaan` fit in the supported scope.
#' @param weight The matrix in the quadratic form, as in [semtest()]:
#'   `"curvature"` for the selected saturated curvature (an observed one must be
#'   positive in the tested directions), `"expected_curvature"` for expected
#'   saturated curvature whatever the curvature used for projection, or
#'   `"score_covariance"` for full efficient-score studentization.
#' @param score_covariance `"uncentered"` for raw score products, or `"centered"`
#'   for separate within-group centering.
#' @param keep Retain normalized matrices, tangent and projection.
#' @param tolerance Matrix validation tolerance.
#' @param curvature Expected or observed saturated curvature for projection.
#' @return A `semTests_test` recording the hypothesis, sensitivity, covariance
#'   convention, group sizes, source version and nuisance stationarity diagnostic.
#' @keywords internal
#' @examples
#' fit <- lavaan::cfa("f =~ x1 + x2 + x3 + x4",
#'                    lavaan::HolzingerSwineford1939, estimator = "ML")
#' s <- semTests:::score_test(fit)
#' pvalues(s, c("SB", "ALL", "PEBA2"))
score_test <- function(object,
                        weight = c("curvature", "expected_curvature", "score_covariance"),
                        score_covariance = c("uncentered", "centered", "expected"), keep = FALSE,
                        tolerance = 1e-10, curvature = c("expected", "observed")) {
  check_ml_inference_adapter(object, "score_test()", allow_fiml = TRUE)
  curvature <- match.arg(curvature)
  weight <- match.arg(weight)
  score_covariance <- match.arg(score_covariance)
  inference_options("Exact fit against saturated moments", keep, tolerance)
  # Only weight = "curvature" with observed curvature weights by the observed
  # matrix; every other choice weights by expected curvature or studentizes.
  observed_weight <- weight == "curvature" && curvature == "observed"
  # The reduced complete-data route has no expected-covariance slot; the pattern
  # route carries the per-pattern expected information needed for it.
  parts <- if (!is_fiml(object) && curvature == "expected" &&
               score_covariance != "expected") {
    ml_global_score_parts(object, score_covariance)
  } else {
    gaussian_global_score_parts(object, score_covariance, curvature, observed_weight)
  }
  result <- score_test_matrices(parts$score,
    metric = if (weight == "score_covariance") NULL else parts$metric,
    score_covariance = parts$score_covariance, nobs = parts$nobs,
    hypothesis = "Exact fit against saturated moments", keep = keep, tolerance = tolerance)
  result$provenance <- c(result$provenance, list(lavaan_version = as.character(utils::packageVersion("lavaan")),
    estimator = requested_estimator(object), curvature = curvature, weight = weight,
    score_covariance = score_covariance, group_nobs = parts$group_nobs))
  result$provenance$backend <- if (is_fiml(object)) "lavaan FIML pattern scores" else "lavaan complete-data ML"
  result$diagnostics$nuisance_gradient <- parts$nuisance_gradient
  result$diagnostics$nuisance_stationarity_norm <- max(abs(parts$nuisance_gradient), 0)
  if (keep) {
    result$ingredients$projection <- parts$projection
    result$ingredients$tangent <- parts$tangent
    result$ingredients$information_blocks <- parts$information
  }
  result
}
