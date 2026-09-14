# Matrix-only inference. Extraction, nuisance profiling, and the scientific
# hypothesis belong to the caller; these kernels do not fit a model.
inference_vector <- function(x, name) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || any(!is.finite(x))) {
    semtests_abort(paste0(name, " must be a nonempty finite numeric vector."),
                  "semTests_error_invalid_inference")
  }
  as.numeric(x)
}

inference_symmetric <- function(x, size, name, tolerance) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(size, size)) ||
      any(!is.finite(x))) {
    semtests_abort(paste0(name, " must be a finite ", size, " by ", size, " matrix."),
                  "semTests_error_invalid_inference")
  }
  scale <- max(abs(x))
  if (max(abs(x - t(x))) > tolerance * scale) {
    semtests_abort(paste0(name, " must be symmetric."), "semTests_error_invalid_inference")
  }
  (x + t(x)) / 2
}

inference_chol <- function(x, name) {
  if (any(!is.finite(x))) {
    semtests_abort(paste0(name, " contains non-finite entries after normalization."),
                  "semTests_error_invalid_inference")
  }
  out <- tryCatch(chol(x), error = function(e) NULL)
  if (is.null(out)) {
    semtests_abort(paste0(name, " must be positive definite; singular or indefinite ",
                         "inference is not supported by this kernel."),
                  "semTests_error_unstable_spectrum")
  }
  out
}

inference_options <- function(hypothesis, keep, tolerance) {
  if (length(hypothesis) != 1L || !is.character(hypothesis) ||
      is.na(hypothesis) || !nzchar(hypothesis) ||
      length(keep) != 1L || !is.logical(keep) || is.na(keep) ||
      length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) ||
      tolerance <= 0 || tolerance >= 1) {
    semtests_abort("Supply a nonempty hypothesis, logical keep, and tolerance in (0, 1).",
                  "semTests_error_invalid_inference")
  }
}

new_inference_test <- function(type, statistic, eigenvalues, reference, hypothesis,
                               normalization, nobs, diagnostics, ingredients) {
  df <- length(eigenvalues)
  p <- if (reference == "chisq") {
    stats::pchisq(statistic, df, lower.tail = FALSE)
  } else {
    pvalue_all(statistic, eigenvalues)
  }
  structure(list(type = type, statistic = statistic, df = df, rank = df,
                 eigenvalues = eigenvalues, reference = reference, p_value = p,
                 hypothesis = hypothesis, normalization = normalization, nobs = nobs,
                 diagnostics = diagnostics, ingredients = ingredients,
                 provenance = list(backend = "supplied matrices",
                                   semTests_version = as.character(utils::packageVersion("semTests")))),
            class = "semTests_test")
}

#' Score test from efficient-score matrices
#'
#' Construct a quadratic score test from an already nuisance-adjusted score,
#' metric, and covariance. No model fitting or nuisance projection is performed.
#' The caller must establish the target hypothesis, correct derivatives, and
#' covariance estimator, including estimated-weight effects where applicable.
#'
#' With summed inputs, the statistic is `score' solve(metric, score)` and the
#' reference weights are the eigenvalues of the covariance whitened by the
#' metric. Set `metric = NULL` to use the covariance itself, giving a fully
#' studentized score test with a chi-square reference. All matrices must be
#' positive definite; this initial interface rejects singular inference.
#'
#' For `normalization = "mean"`, supply the average score and per-observation
#' metric and score covariance. The covariance is **not** the covariance of the
#' average score. `nobs` scales all three inputs to summed units. Centering the
#' score covariance is the caller's choice and must match the intended method.
#'
#' [pvalues()] applies existing eigenvalue calibrations to a supplied metric
#' score test. Its default here is `"ALL"`, the unpenalized mixture law.
#' Penalized calibrations are available for research comparisons; numerical
#' agreement does not establish their calibration for a new score procedure.
#' Fully studentized score and Wald results accept only `"STD"`.
#'
#' @param score Nonempty numeric efficient-score vector in tested coordinates.
#' @param metric Positive-definite quadratic metric, or `NULL` for full
#'   covariance studentization.
#' @param score_covariance Covariance of the summed score (summed inputs), or
#'   per-observation score covariance (mean inputs).
#' @param normalization Either `"sum"` or `"mean"`.
#' @param nobs Total sample size. Required for mean inputs, optional for sums.
#' @param hypothesis Description of the hypothesis represented by the inputs.
#' @param keep Keep the normalized input matrices in the result for inspection.
#' @param tolerance Relative tolerance for matrix symmetry and PSD validation.
#' @return A `semTests_test` list containing statistic, df, rank, eigenvalues,
#'   reference law, default p-value, hypothesis, normalization, diagnostics,
#'   provenance and optional ingredients. Degrees of freedom are the dimension
#'   of the supplied tested coordinates. No misspecification-validity claim is
#'   inferred from the supplied matrices.
#' @keywords internal
#' @examples
#' s <- semTests:::score_test_matrices(c(1, 2), diag(2), diag(c(1, 2)))
#' pvalues(s, c("SB", "ALL", "PEBA2"))
#' as.data.frame(s)
score_test_matrices <- function(score, metric = NULL, score_covariance,
                                 normalization = c("sum", "mean"), nobs = NULL,
                                 hypothesis = "Supplied efficient-score hypothesis",
                                 keep = FALSE, tolerance = 1e-10) {
  inference_options(hypothesis, keep, tolerance)
  normalization <- match.arg(normalization)
  score <- inference_vector(score, "score")
  d <- length(score)
  score_covariance <- inference_symmetric(score_covariance, d, "score_covariance", tolerance)
  if (!is.null(nobs) && (length(nobs) != 1L || !is.numeric(nobs) ||
      !is.finite(nobs) || nobs < 1 || nobs != floor(nobs))) {
    semtests_abort("nobs must be a positive integer.", "semTests_error_invalid_inference")
  }
  if (normalization == "mean" && is.null(nobs)) {
    semtests_abort("Mean inputs require nobs.", "semTests_error_invalid_inference")
  }
  studentized <- is.null(metric)
  metric <- if (studentized) score_covariance else inference_symmetric(metric, d, "metric", tolerance)
  if (normalization == "mean") {
    score <- nobs * score
    metric <- nobs * metric
    score_covariance <- nobs * score_covariance
  }
  # Both matrices need full rank, including for the mixture comparator.
  inference_chol(score_covariance, "score_covariance")
  R <- inference_chol(metric, "metric")
  z <- forwardsolve(t(R), score)
  statistic <- sum(z^2)
  if (!is.finite(statistic)) {
    semtests_abort("Score quadratic overflowed.", "semTests_error_invalid_inference")
  }
  eigenvalues <- if (studentized) rep(1, d) else {
    left <- forwardsolve(t(R), score_covariance)
    K <- t(forwardsolve(t(R), t(left)))
    values <- eigen((K + t(K)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(values)) || any(values <= 0)) {
      semtests_abort("Whitened covariance is numerically singular or indefinite.",
                    "semTests_error_unstable_spectrum")
    }
    values
  }
  new_inference_test("score", statistic, eigenvalues,
    if (studentized) "chisq" else "weighted_chisq", hypothesis,
    normalization = normalization, nobs = nobs,
    diagnostics = list(metric_rcond = rcond(metric),
                       score_covariance_rcond = rcond(score_covariance),
                       tolerance = tolerance, studentized = studentized),
    ingredients = if (keep) list(score = score, metric = metric, score_covariance = score_covariance,
                                 units = "sum") else NULL)
}

#' Linear Wald test from parameter estimates and sampling covariance
#'
#' Test `R %*% theta = rhs` using the supplied **sampling covariance** of theta.
#' The statistic is the restriction contrast whitened by `R %*% vcov %*% t(R)`.
#' There is no additional sample-size multiplication. To use an asymptotic
#' covariance for `sqrt(n) * theta`, divide it by `n` before calling this function.
#'
#' The full covariance may be positive semidefinite, as with constrained
#' parameterizations, but the restriction covariance must be positive definite.
#' Redundant restrictions are rejected. Inference under fixed misspecification
#' requires the caller to supply the complete sandwich covariance for the actual
#' pseudo-parameter. This function does not estimate it.
#'
#' @param theta Numeric vector of parameter estimates.
#' @param R Numeric restriction matrix, one row per tested restriction, with
#'   columns in the same positional order as theta and vcov.
#' @param vcov Sampling covariance matrix in theta's positional order.
#' @param rhs Numeric vector of null values, one per restriction; `NULL` means zeros.
#' @inheritParams score_test_matrices
#' @return A `semTests_test` result with a chi-square reference; see
#'   [score_test_matrices()] for its fields.
#' @keywords internal
#' @examples
#' w <- semTests:::wald_test_matrices(c(.4, .2), matrix(c(1, -1), 1), diag(c(.02, .03)))
#' pvalues(w)
wald_test_matrices <- function(theta, R, vcov, rhs = NULL,
                                hypothesis = "R theta = rhs", keep = FALSE,
                                tolerance = 1e-10) {
  inference_options(hypothesis, keep, tolerance)
  theta <- inference_vector(theta, "theta")
  V <- inference_symmetric(vcov, length(theta), "vcov", tolerance)
  if (!is.matrix(R) || !is.numeric(R) || ncol(R) != length(theta) ||
      !nrow(R) || any(!is.finite(R))) {
    semtests_abort("R must be a finite matrix with one column per parameter.",
                  "semTests_error_invalid_inference")
  }
  row_scale <- sqrt(rowSums(R^2))
  if (any(row_scale == 0) ||
      qr(t(sweep(R, 1L, row_scale, "/")), tol = tolerance)$rank != nrow(R)) {
    semtests_abort("R contains redundant restrictions.", "semTests_error_unstable_spectrum")
  }
  rhs <- if (is.null(rhs)) rep(0, nrow(R)) else inference_vector(rhs, "rhs")
  if (length(rhs) != nrow(R)) {
    semtests_abort("rhs must have one value per restriction.", "semTests_error_invalid_inference")
  }
  values <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (any(values < -tolerance * max(abs(values)))) {
    semtests_abort("vcov must be positive semidefinite.", "semTests_error_invalid_inference")
  }
  contrast <- drop(R %*% theta) - rhs
  restriction_covariance <- R %*% V %*% t(R)
  root <- inference_chol((restriction_covariance + t(restriction_covariance)) / 2,
                         "Restriction covariance")
  statistic <- sum(forwardsolve(t(root), contrast)^2)
  if (!is.finite(statistic)) {
    semtests_abort("Wald quadratic overflowed.", "semTests_error_invalid_inference")
  }
  new_inference_test("wald", statistic, rep(1, nrow(R)), "chisq", hypothesis,
    normalization = "sampling covariance", nobs = NULL,
    diagnostics = list(restriction_rcond = rcond(restriction_covariance), tolerance = tolerance),
    ingredients = if (keep) list(theta = theta, R = R, rhs = rhs, vcov = V,
                                 contrast = contrast, restriction_covariance = restriction_covariance) else NULL)
}

#' @export
print.semTests_test <- function(x, ...) {
  cat(x$type, "test:", x$hypothesis, "\n")
  print(as.data.frame(x), row.names = FALSE, ...)
  invisible(x)
}

#' @export
as.data.frame.semTests_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(type = x$type, statistic = x$statistic, df = x$df, rank = x$rank,
             reference = x$reference, p_value = x$p_value, hypothesis = x$hypothesis,
             row.names = row.names, stringsAsFactors = FALSE)
}

