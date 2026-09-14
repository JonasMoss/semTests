# Computation settings for the inference constructors. Keeping these out of the
# public signatures leaves only the question being asked in `semtest()`, and
# gives one documented place for every convention choice.

#' Computation settings for `semtest()`
#'
#' Collects the convention choices used when building a test. These are
#' deliberately separated from `semtest()`'s own arguments: the signature states
#' the question, the control states how the matrices are formed.
#'
#' `curvature` and `score_covariance` default to `NULL`, meaning each test picks
#' the convention appropriate to the null it evaluates. Goodness of fit assumes
#' the model holds, since that is its hypothesis, and uses expected curvature
#' with uncentered score products. Nested restriction tests allow the baseline
#' model to be wrong and use observed curvature with centered score products.
#' Supplying a value overrides the choice everywhere it applies.
#'
#' The resolved settings are recorded in the returned test's provenance, so a
#' result always reports the conventions that produced it.
#'
#' @param curvature `"observed"` or `"expected"` sensitivity. `NULL` defers to
#'   the test. Expected and observed agree in probability when the model is
#'   correctly specified, and differ otherwise.
#' @param score_covariance `"centered"`, `"uncentered"` or `"expected"` score
#'   covariance. `NULL` defers to the test. Centered subtracts group score
#'   means and stays valid when the score has nonzero mean; uncentered is
#'   consistent only at exact fit. `"expected"` uses model-based information in
#'   the covariance slot; see Details.
#' @param weight_effect `"estimated"` propagates first-stage weight estimation
#'   into the covariance; `"fixed"` treats the weights as known. Ignored for ML.
#' @param derivative_step Relative central-difference step for numerical
#'   criterion curvature. Must lie in (0, .01].
#' @param tolerance Relative tolerance for matrix symmetry and definiteness
#'   validation. Must lie in (0, 1).
#' @param keep Retain the input matrices and projections in the result.
#'
#' @details
#' `score_covariance = "expected"` places model-based information in the
#' covariance slot rather than an empirical score product. Under complete data
#' the information equality makes this the same matrix as the expected
#' curvature, so the weighted spectrum is flat and the test reduces to the
#' uncorrected chi-square already available as `tests = "STD"`. It is a
#' genuinely distinct test only where the information equality fails, which is
#' the case under missing data, where the expectation runs over the missingness
#' patterns as well. That is why it is worth requesting for FIML and not
#' otherwise.
#'
#' @return An object of class `semTests_control`.
#' @export
#' @examples
#' inference_control(curvature = "observed")
#' inference_control(score_covariance = "expected", keep = TRUE)
inference_control <- function(curvature = NULL,
                              score_covariance = NULL,
                              weight_effect = c("estimated", "fixed"),
                              derivative_step = 1e-4,
                              tolerance = 1e-10,
                              keep = FALSE) {
  curvature <- control_level(curvature, c("observed", "expected"), "curvature")
  score_covariance <- control_level(score_covariance,
    c("centered", "uncentered", "expected"), "score_covariance")
  weight_effect <- match.arg(weight_effect)
  if (length(derivative_step) != 1L || !is.numeric(derivative_step) ||
      !is.finite(derivative_step) || derivative_step <= 0 || derivative_step > .01) {
    semtests_abort("derivative_step must be a single number in (0, .01].",
                   "semTests_error_invalid_inference")
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) ||
      tolerance <= 0 || tolerance >= 1) {
    semtests_abort("tolerance must be a single number in (0, 1).",
                   "semTests_error_invalid_inference")
  }
  if (length(keep) != 1L || !is.logical(keep) || is.na(keep)) {
    semtests_abort("keep must be TRUE or FALSE.", "semTests_error_invalid_inference")
  }
  structure(list(curvature = curvature, score_covariance = score_covariance,
                 weight_effect = weight_effect, derivative_step = derivative_step,
                 tolerance = tolerance, keep = keep),
            class = "semTests_control")
}

# NULL means "the test decides"; anything else must be one of the levels. This
# is not match.arg, which would silently accept the full default vector.
control_level <- function(value, levels, name) {
  if (is.null(value)) return(NULL)
  if (length(value) != 1L || !is.character(value) || is.na(value) ||
      !value %in% levels) {
    semtests_abort(sprintf("%s must be NULL or one of %s.", name,
                           paste(sprintf('"%s"', levels), collapse = ", ")),
                   "semTests_error_invalid_inference")
  }
  value
}

#' @export
print.semTests_control <- function(x, ...) {
  shown <- vapply(x, function(v) if (is.null(v)) "<by test>" else as.character(v),
                  character(1))
  cat("semTests inference control\n")
  cat(paste0("  ", format(names(shown)), "  ", shown, collapse = "\n"), "\n", sep = "")
  invisible(x)
}

# Fill unset conventions from the null being tested. Goodness of fit assumes the
# model holds; restriction tests do not. Keep the two default sets here so the
# divergence is visible in one place rather than spread over signatures.
resolve_control <- function(control, context) {
  if (!inherits(control, "semTests_control")) {
    semtests_abort("control must come from inference_control().",
                   "semTests_error_invalid_inference")
  }
  defaults <- if (identical(context, "fit")) {
    list(curvature = "expected", score_covariance = "uncentered")
  } else {
    list(curvature = "observed", score_covariance = "centered")
  }
  if (is.null(control$curvature)) control$curvature <- defaults$curvature
  if (is.null(control$score_covariance)) {
    control$score_covariance <- defaults$score_covariance
  }
  control
}

# Merge `...` into a control, rejecting names that are not control settings.
# Silently dropping a typo here would silently change the conventions.
control_dots <- function(control, dots) {
  if (!length(dots)) return(control)
  known <- names(formals(inference_control))
  unknown <- setdiff(names(dots), known)
  if (length(unknown) || is.null(names(dots)) || any(!nzchar(names(dots)))) {
    semtests_abort(sprintf(
      "Unknown inference setting%s: %s. Settings are: %s.",
      if (length(unknown) == 1L) "" else "s",
      paste(unknown, collapse = ", "), paste(known, collapse = ", ")),
      "semTests_error_invalid_inference")
  }
  do.call(inference_control, utils::modifyList(as.list(control)[known], dots))
}
