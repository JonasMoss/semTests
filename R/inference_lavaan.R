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
#' misspecification, use [score_test_nested()] or [wald_test()].
#'
#' Supports normal-likelihood ML/MLM/MLR and FIML, one or several groups,
#' affine equalities, unweighted raw data, and joint random-x inference.
#' Categorical data, sampling weights, clusters, bounds and nonlinear
#' constraints are excluded. A fully studentized score uses a chi-square
#' reference. Other metrics use the full eigenvalue mixture by default.
#'
#' @param object A converged `lavaan` fit in the supported scope.
#' @param metric `"information"` for the Fisher-metric quadratic, or
#'   `"sandwich"` for full efficient-score studentization, or `"observed"`
#'   for observed saturated curvature (must be positive in tested directions).
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
score_test <- function(object, metric = c("information", "sandwich", "observed"),
                        score_covariance = c("uncentered", "centered", "expected"), keep = FALSE,
                        tolerance = 1e-10, curvature = c("expected", "observed")) {
  check_ml_inference_adapter(object, "score_test()", allow_fiml = TRUE)
  curvature <- match.arg(curvature)
  metric <- match.arg(metric)
  score_covariance <- match.arg(score_covariance)
  inference_options("Exact fit against saturated moments", keep, tolerance)
  # The reduced complete-data route has no expected-covariance slot; the pattern
  # route carries the per-pattern expected information needed for it.
  parts <- if (!is_fiml(object) && curvature == "expected" && metric != "observed" &&
               score_covariance != "expected") {
    ml_global_score_parts(object, score_covariance)
  } else {
    gaussian_global_score_parts(object, score_covariance, curvature, metric)
  }
  result <- score_test_matrices(parts$score,
    metric = if (metric == "sandwich") NULL else parts$metric,
    score_covariance = parts$score_covariance, nobs = parts$nobs,
    hypothesis = "Exact fit against saturated moments", keep = keep, tolerance = tolerance)
  result$provenance <- c(result$provenance, list(lavaan_version = as.character(utils::packageVersion("lavaan")),
    estimator = requested_estimator(object), curvature = curvature, metric = metric,
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

#' Linear Wald restrictions with fitted or complete-sandwich covariance
#'
#' Test `R theta = rhs`. The default `vcov = "fitted"` uses lavaan's stored
#' sampling covariance for ML/FIML, preserving its SE/information conventions.
#' `vcov = "sandwich"` constructs a new criterion-based sandwich. Its default
#' observed sensitivity includes full model curvature. Gaussian case scores
#' include fitted-mean residual contributions; LS includes joint moment and
#' estimated-weight influences, with their cross covariance.
#'
#' Complete-sandwich inference supports ML/FIML, continuous GLS/WLS/DWLS/ULS,
#' and complete all-ordered DWLS/WLS/ULS, with affine equalities and unweighted
#' raw observations. The ordinal path differentiates lavaan's native empirical
#' frequency functional, refines polychoric roots, and uses Richardson
#' extrapolation. It is a numerical reference implementation, potentially slow
#' with many distinct response patterns. Empty pairwise cells, mixed indicators,
#' sampling weights, clusters and nonlinear constraints are excluded.
#'
#' The ordinal weight map retains lavaan's score-product curvature. The influence
#' of its first-stage estimator uses actual derivative sensitivity. These are
#' distinct under first-stage misspecification; see the developer validation
#' notes for the corresponding magmaan comparison.
#'
#' Fixed misspecification inference concerns the criterion's pseudo-parameter,
#' which need not equal a generating-model parameter. `curvature = "expected"`
#' and `weight_effect = "fixed"` are explicit comparators, not complete
#' corrections. There is no second eigenvalue correction after studentization.
#'
#' @param object A converged supported `lavaan` fit.
#' @param R Restriction matrix with columns in positional `lavaan::coef(object)`
#'   order. Names are not used to reorder columns.
#' @param rhs Null restriction values, or `NULL` for zero.
#' @param hypothesis Description of the linear hypothesis and population target.
#' @param keep Retain parameter covariance and estimating-equation matrices.
#' @param tolerance Matrix validation tolerance.
#' @param vcov `"fitted"` or a newly constructed `"sandwich"` covariance.
#' @param curvature Observed or expected criterion sensitivity for the new
#'   sandwich. Ignored when `vcov = "fitted"`.
#' @param score_covariance Group-centered or uncentered Gaussian score products.
#'   LS requires `"centered"`. Ignored when `vcov = "fitted"`.
#' @param weight_effect Include estimated weights or treat them as fixed when
#'   constructing an LS sandwich. Ignored for ML and fitted covariance.
#' @param derivative_step Relative central-difference step for criterion
#'   curvature. The ordinal probability perturbation uses its square root times
#'   the response-pattern frequency, with Richardson extrapolation.
#' @return A `semTests_test` with chi-square reference, parameter order and
#'   covariance conventions. Complete-sandwich results also report stationarity
#'   and numerical curvature asymmetry. No parameter refitting is performed.
#' @keywords internal
#' @examples
#' fit <- lavaan::cfa("f =~ x1 + a*x2 + b*x3 + x4",
#'                    lavaan::HolzingerSwineford1939, estimator = "MLM")
#' R <- matrix(0, 1, length(lavaan::coef(fit)))
#' R[1, match(c("a", "b"), names(lavaan::coef(fit)))] <- c(1, -1)
#' pvalues(semTests:::wald_test(fit, R))
wald_test <- function(object, R, rhs = NULL, hypothesis = "R theta = rhs",
                       keep = FALSE, tolerance = 1e-10,
                       vcov = c("fitted", "sandwich"),
                       curvature = c("observed", "expected"),
                       score_covariance = c("centered", "uncentered", "expected"),
                       weight_effect = c("estimated", "fixed"),
                       derivative_step = 1e-4) {
  vcov <- match.arg(vcov)
  curvature <- match.arg(curvature)
  score_covariance <- match.arg(score_covariance)
  weight_effect <- match.arg(weight_effect)
  if (vcov == "fitted") check_ml_inference_adapter(object, "wald_test()", allow_fiml = TRUE)
  else check_criterion_adapter(object)
  if (vcov == "sandwich" && object@Options$estimator != "ML" && score_covariance != "centered")
    stop("LS requires centered first-stage influence functions.", call. = FALSE)
  inference_options(hypothesis, keep, tolerance)
  theta <- lavaan::coef(object)
  if (is.matrix(R) && is.numeric(R) && ncol(R) == length(theta) &&
      nrow(R) > 0L && all(is.finite(R))) {
    scales <- sqrt(rowSums(R^2))
    if (all(scales > 0)) {
      reduced <- sweep(R, 1L, scales, "/") %*% model_parameter_basis(object)
      singular <- svd(reduced, nu = 0, nv = 0)$d
      if (sum(singular > tolerance * max(1, singular)) != nrow(R)) {
        semtests_abort("R contains restrictions already fixed or dependent in the fitted model.",
                      "semTests_error_unstable_spectrum")
      }
    }
  }
  parts <- NULL
  V <- if (vcov == "fitted") {
    tryCatch(lavaan::lavInspect(object, "vcov"), error = function(e) NULL)
  } else {
    parts <- parameter_inference_parts(object, curvature, score_covariance, weight_effect, derivative_step)
    inverse <- inference_solve(parts$curvature)
    B <- model_parameter_basis(object)
    B %*% inverse %*% parts$score_covariance %*% inverse %*% t(B)
  }
  if (is.null(V) || !length(V)) {
    semtests_abort("wald_test() requires a fitted sampling covariance; refit with standard errors.",
                  "semTests_error_invalid_inference")
  }
  result <- wald_test_matrices(theta, R, as.matrix(V), rhs,
    hypothesis = hypothesis, keep = keep, tolerance = tolerance)
  result$nobs <- object@SampleStats@ntotal
  result$provenance <- c(result$provenance, list(lavaan_version = as.character(utils::packageVersion("lavaan")),
    estimator = requested_estimator(object), se = object@Options$se,
    information = object@Options$information, parameter_names = names(theta)))
  result$provenance$backend <- if (vcov == "fitted") "lavaan fitted covariance" else "semTests complete criterion sandwich"
  if (!is.null(parts)) {
    result$provenance <- c(result$provenance, parts$provenance)
    result$diagnostics$stationarity_norm <- max(abs(parts$score))/parts$nobs
    result$diagnostics$curvature_asymmetry <- parts$curvature_asymmetry
    if (keep) result$ingredients$estimating_equations <- parts
  }
  result
}
