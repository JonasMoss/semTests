# Defensive checks on the fitted object, kept separate from the p-value engine.

#' Warn when a FIML fit uses expected information.
#'
#' Under missing data the expected (Fisher) information is consistent only under
#' MCAR; under MAR -- the regime FIML is adopted for -- the observed information
#' is required for valid inference (Kenward & Molenberghs, 1998). lavaan defaults
#' to observed information for `missing = "ml"`, but it silently accepts
#' `information = "expected"`, which would make the eigenvalue p-values rest on
#' the stronger MCAR assumption. This emits a warning (not an error: expected is
#' still valid under MCAR).
#'
#' lavaan stores `information` as a length-2 vector, so we test its first entry.
#'
#' @param fit A fitted `lavaan` object.
#' @return `fit`, invisibly.
#' @keywords internal
warn_fiml_information <- function(fit) {
  if (!identical(fit@Options$missing, "listwise") &&
      identical(fit@Options$information[1], "expected")) {
    warning(
      "Missing-data (FIML) fit with expected information, which is valid only ",
      "under MCAR. Refit with information = \"observed\" (lavaan's default for ",
      "FIML) for inference valid under MAR.",
      call. = FALSE
    )
  }
  invisible(fit)
}

# Estimators whose minimum-discrepancy statistic the eigenvalue correction is
# valid for. lavaan reports the *base* estimator here, so the MLM/MLR/WLSMV
# shortcuts all collapse to ML/DWLS; the supported surface is documented in
# `?semTests-support`.
.supported_estimators <- c("ML", "GLS", "ULS", "DWLS", "WLS")

#' Reject a fit whose configuration is outside the supported surface.
#'
#' The single entry-level gate for [pvalues()]. It admits exactly the
#' configurations documented in [semTests-support] -- the supported continuous
#' and categorical estimators, complete data, and single-group FIML -- and stops
#' with a pointer to `?semTests-support` otherwise. Statistic- and gamma-specific
#' rejections (the normal-theory-only RLS statistic and `UG` gamma) depend on the
#' parsed test string and stay with the code that consumes them (`make_chisqs()`,
#' [gamma()]); this gate is purely fit-shape.
#'
#' @param fit A fitted `lavaan` object.
#' @return `fit`, invisibly.
#' @keywords internal
check_supported <- function(fit) {
  est <- fit@Options$estimator
  if (!est %in% .supported_estimators) {
    stop("Estimator '", est, "' is not supported. The supported estimators are ",
         "ML/MLM/MLR, GLS, ULS, and categorical WLSMV/DWLS (ADF/WLS is a ",
         "degenerate exception). See `?semTests-support`.", call. = FALSE)
  }
  if (is_fiml(fit)) {
    # Continuous by construction (is_fiml excludes categorical); enforce the
    # single-group / no-fixed-exogenous-covariate FIML constraints here so they
    # surface at the entry point rather than deep in the FIML engine.
    fiml_check_supported(fit, "FIML (missing data)")
  } else if (!identical(fit@Options$missing, "listwise")) {
    stop("Missing data is only supported through FIML (continuous ML/MLR with ",
         "missing = \"ml\"); this fit uses missing = \"", fit@Options$missing[1],
         "\". See `?semTests-support`.", call. = FALSE)
  }
  invisible(fit)
}

#' Reject a nested pair whose configuration is outside the supported surface.
#'
#' The single entry-level gate for [pvalues_nested()]: categorical nesting is
#' deferred, missing-data nesting requires both fits to be FIML, and FIML nesting
#' supports `method = "2000"` only. Each fit is also run through
#' [check_supported()]. The `UG`-gamma rejection for FIML stays in the p-value
#' engine because it depends on the parsed test string.
#'
#' @param m0,m1 Two nested `lavaan` objects (canonical order: `m0` constrained).
#' @param method The nested reduction method, `"2000"` or `"2001"`.
#' @param A.method The FIML restriction map, `"exact"` or `"delta"` (validated
#'   upstream; accepted here for signature completeness).
#' @return `TRUE`, invisibly.
#' @keywords internal
check_supported_nested <- function(m0, m1, method, A.method = "exact") {
  if (isTRUE(m0@Model@categorical) || isTRUE(m1@Model@categorical)) {
    stop("Nested tests are not yet supported for categorical data; ",
         "single-model tests are. See `?semTests-support`.", call. = FALSE)
  }
  check_supported(m0)
  check_supported(m1)
  if ((!identical(m0@Options$missing, "listwise") ||
       !identical(m1@Options$missing, "listwise")) &&
      !(is_fiml(m0) && is_fiml(m1))) {
    stop("Nested tests for missing data currently require FIML fits. ",
         "See `?semTests-support`.", call. = FALSE)
  }
  # If either fit is FIML then both are: the mixed-missing check above rejects a
  # FIML/non-FIML pair, and check_supported() rejects any other missing mode. So
  # only the method restriction remains to enforce here.
  if (is_fiml(m0) || is_fiml(m1)) {
    if (method != "2000") {
      stop("Nested FIML tests support method = \"2000\" only. ",
           "See `?semTests-support`.", call. = FALSE)
    }
  }
  invisible(TRUE)
}
