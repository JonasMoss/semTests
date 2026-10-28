# Single public constructor for the score and criterion tests. The internal
# adapters keep their own signatures; this states the question and resolves
# conventions once.

#' Build a score or criterion test
#'
#' Constructs an inspectable test object for overall fit or for a nested
#' restriction. Pass the result to [pvalues()] for calibrated p-values.
#'
#' Supply `alternative` to test a restriction, or leave it `NULL` to test
#' overall fit against saturated means and covariances. The restricted model
#' comes first, as in [pvalues_nested()]. The alternative may be created with
#' `do.fit = FALSE` when `statistic = "score"`: it supplies the
#' parameterization and is evaluated at the restricted estimates without
#' refitting.
#'
#' `statistic = "score"` evaluates the score at the restricted fit and needs
#' only that fit to be optimized. `statistic = "criterion"` uses the actual
#' difference of the summed fitting criteria and needs both fits. Their null is
#' a common population minimizer under the shared criterion, which allows the
#' baseline model to be misspecified. That is a different null from
#' [pvalues_nested()], which tests exact fit of both models using Satorra's
#' 2000 construction, so the two need not agree.
#'
#' Overall fit with `statistic = "criterion"` is the classical discrepancy
#' statistic, which is what [pvalues()] computes; `semtest()` does not duplicate
#' it.
#'
#' @param object A converged `lavaan` fit. For a restriction test this is the
#'   restricted model.
#' @param alternative A nested alternative fit or unoptimized template, using
#'   the same observations, estimator and moment conventions. `NULL` tests
#'   overall fit.
#' @param statistic `"score"` or `"criterion"`. `"criterion"` requires
#'   `alternative`.
#' @param weight Which matrix weights the quadratic form, and therefore which
#'   statistic is computed. This is a separate decision from `curvature`, which
#'   selects the sensitivity used to project out nuisance directions and so
#'   enters the reference law.
#'
#'   * `"curvature"` weights by whichever sensitivity `curvature` selected.
#'   * `"expected_curvature"` weights by expected curvature even when
#'     `curvature = "observed"`. This is the combination that keeps observed
#'     sensitivity in the projection, which is what the limit requires when the
#'     model may be misspecified, while weighting by the more stable expected
#'     matrix.
#'   * `"score_covariance"` studentizes fully and gives a chi-square reference.
#'
#'   The first two give a weighted chi-square reference, so the calibrations in
#'   [pvalues()] apply. Ignored when `statistic = "criterion"`, which weights by
#'   curvature by construction.
#' @param control Conventions from [inference_control()].
#' @param ... Individual [inference_control()] settings, merged into `control`.
#'   Unknown names are an error.
#' @return A `semTests_test`, for use with [pvalues()], `print()` and
#'   `as.data.frame()`.
#' @seealso [pvalues()] for calibrated p-values and for the classical
#'   discrepancy statistics, [inference_control()] for the conventions.
#' @export
#' @examples
#' model <- "f =~ x1 + a*x2 + b*x3 + x4"
#' d <- lavaan::HolzingerSwineford1939
#' fit <- lavaan::cfa(model, d, meanstructure = TRUE)
#'
#' # Overall fit, score statistic.
#' pvalues(semtest(fit), c("SB", "ALL"))
#'
#' # Can the equality be released?
#' restricted <- lavaan::cfa(paste(model, "a == b", sep = "\n"), d,
#'                           meanstructure = TRUE)
#' semtest(restricted, fit, statistic = "score", weight = "score_covariance")
#' semtest(restricted, fit, statistic = "criterion")
semtest <- function(object, alternative = NULL,
                    statistic = c("score", "criterion"),
                    weight = c("curvature", "expected_curvature", "score_covariance"),
                    control = inference_control(), ...) {
  statistic <- match.arg(statistic)
  weight <- match.arg(weight)
  control <- resolve_control(control_dots(control, list(...)),
                             if (is.null(alternative)) "fit" else "nested")
  if (is.null(alternative)) {
    if (statistic == "criterion") {
      semtests_abort(paste0("Overall fit with statistic = \"criterion\" is the ",
        "classical discrepancy statistic; use pvalues() for it."),
        "semTests_error_invalid_inference")
    }
    return(semtest_fit(object, weight, control))
  }
  if (statistic == "criterion") {
    # An unfitted template carries a parameterization but no criterion minimum.
    # Say that, rather than reporting it as an estimation failure.
    if (inherits(alternative, "lavaan") && !alternative@optim$converged &&
        identical(alternative@optim$iterations, 0L)) {
      semtests_abort(paste0(
        "statistic = \"criterion\" needs the alternative's fitting criterion, ",
        "so it must be optimized. An unoptimized template works with ",
        "statistic = \"score\"."), "semTests_error_invalid_inference")
    }
    return(criterion_test_nested(object, alternative,
      curvature = control$curvature, score_covariance = control$score_covariance,
      weight_effect = control$weight_effect,
      derivative_step = control$derivative_step,
      keep = control$keep, tolerance = control$tolerance))
  }
  score_test_nested(object, alternative, weight = weight,
    curvature = control$curvature, score_covariance = control$score_covariance,
    weight_effect = control$weight_effect,
    derivative_step = control$derivative_step,
    keep = control$keep, tolerance = control$tolerance)
}

semtest_fit <- function(object, weight, control) {
  score_test(object, weight = weight,
             score_covariance = control$score_covariance,
             keep = control$keep, tolerance = control$tolerance,
             curvature = control$curvature)
}
