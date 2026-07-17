# Defensive checks on the fitted object, kept separate from the p-value engine.

#' Reject anything that is not a fitted lavaan object.
#'
#' The class gate for [pvalues()] / [pvalues_nested()]. It must run before any
#' `@`-slot access so a `NULL` / data.frame / list argument fails with a readable
#' message instead of a cryptic S4 "no slot of name ..." error -- the nested
#' entry point in particular reads `m0@test` to compute the degrees of freedom
#' before the support gate runs.
#'
#' @param x The object passed by the user.
#' @param arg The argument name, used in the message (e.g. `"object"`, `"m0"`).
#' @return `x`, invisibly.
#' @keywords internal
check_lavaan <- function(x, arg = "object") {
  if (!inherits(x, "lavaan")) {
    what <- if (is.null(x)) {
      "NULL"
    } else {
      paste0("an object of class ",
             paste(dQuote(class(x), q = FALSE), collapse = "/"))
    }
    stop("`", arg, "` must be a fitted lavaan object (class \"lavaan\"), not ",
         what, ". Fit your model with lavaan::cfa()/sem()/lavaan() first. ",
         "See `?semTests-support`.", call. = FALSE)
  }
  invisible(x)
}

# Estimators whose minimum-discrepancy statistic the eigenvalue correction is
# valid for. lavaan reports the *base* estimator here, so the MLM/MLR/WLSMV
# shortcuts all collapse to ML/DWLS; the supported surface is documented in
# `?semTests-support`.
.supported_estimators <- c("ML", "GLS", "ULS", "DWLS", "WLS")
.supported_categorical_estimators <- c("DWLS", "ULS", "WLS")

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
  check_lavaan(fit)
  est <- fit@Options$estimator
  if (!est %in% .supported_estimators) {
    stop("Estimator '", est, "' is not supported. The supported estimators are ",
         "ML/MLM/MLR, GLS, ULS, and categorical DWLS/ULS/WLS families ",
         "(ADF/WLS is a ",
         "degenerate exception). See `?semTests-support`.", call. = FALSE)
  }
  if (isTRUE(fit@Model@categorical)) {
    if (!est %in% .supported_categorical_estimators) {
      stop("Categorical fits currently support lavaan's DWLS, ULS, and WLS ",
           "estimator families; this fit uses '", est,
           "'. See `?semTests-support`.", call. = FALSE)
    }
    if (!fit@Options$missing[1] %in% c("listwise", "pairwise")) {
      stop("Categorical fits support missing = \"listwise\" or \"pairwise\"; ",
           "this fit uses missing = \"", fit@Options$missing[1],
           "\". See `?semTests-support`.", call. = FALSE)
    }
  } else if (is_fiml(fit)) {
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

#' Verify that two categorical fits define one comparable nested problem.
#' @keywords internal
check_categorical_nested_pair <- function(m0, m1) {
  requested_estimator <- function(fit) {
    out <- fit@Options$estimator.orig
    if (is.null(out) || !length(out) || is.na(out[1])) {
      fit@Options$estimator[1]
    } else {
      out[1]
    }
  }
  same <- function(x, y) {
    isTRUE(all.equal(x, y, tolerance = 0, check.attributes = TRUE))
  }

  if (!identical(requested_estimator(m0), requested_estimator(m1))) {
    stop("Nested categorical fits must use the same requested lavaan ",
         "estimator.", call. = FALSE)
  }
  if (!identical(m0@Options$parameterization[1],
                 m1@Options$parameterization[1])) {
    stop("Nested categorical fits must use the same parameterization.",
         call. = FALSE)
  }
  if (!identical(m0@Options$missing[1], m1@Options$missing[1])) {
    stop("Nested categorical fits must use the same missing-data mode.",
         call. = FALSE)
  }
  if (!identical(lavaan::lavNames(m0, "ov"),
                 lavaan::lavNames(m1, "ov")) ||
      !identical(lavaan::lavNames(m0, "ov.ord"),
                 lavaan::lavNames(m1, "ov.ord"))) {
    stop("Nested categorical fits must use the same observed and ordered ",
         "variables in the same order.", call. = FALSE)
  }
  if (!identical(m0@Data@group.label, m1@Data@group.label) ||
      !identical(as.integer(m0@SampleStats@nobs),
                 as.integer(m1@SampleStats@nobs))) {
    stop("Nested categorical fits must use the same groups and group sample ",
         "sizes.", call. = FALSE)
  }
  if (!same(lavaan::lavInspect(m0, "data"),
            lavaan::lavInspect(m1, "data"))) {
    stop("Nested categorical fits must use the same raw data and missingness ",
         "mask.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Reject a nested pair whose configuration is outside the supported surface.
#'
#' The single entry-level gate for [pvalues_nested()]: categorical nesting uses
#' Satorra 2000 with a delta restriction map, missing-data nesting otherwise
#' requires both fits to be FIML, and FIML nesting supports `method = "2000"`
#' only. Each fit is also run through
#' [check_supported()]. The `UG`-gamma rejection for FIML stays in the p-value
#' engine because it depends on the parsed test string.
#'
#' @param m0,m1 Two nested `lavaan` objects (canonical order: `m0` constrained).
#' @param method The nested reduction method, `"2000"` or `"2001"`.
#' @param A.method The FIML restriction map, `"exact"` or `"delta"` (validated
#'   upstream; accepted here for signature completeness).
#' @return `TRUE`, invisibly.
#' @keywords internal
check_supported_nested <- function(m0, m1, method, A.method = "delta") {
  categorical0 <- isTRUE(m0@Model@categorical)
  categorical1 <- isTRUE(m1@Model@categorical)
  if (xor(categorical0, categorical1)) {
    stop("Nested tests require both fits to use the same data type; one fit is ",
         "categorical and the other is continuous. See `?semTests-support`.",
         call. = FALSE)
  }
  check_supported(m0)
  check_supported(m1)
  if (categorical0 && categorical1) {
    if (method != "2000") {
      stop("Nested categorical tests support method = \"2000\" only. ",
           "See `?semTests-support`.", call. = FALSE)
    }
    if (A.method != "delta") {
      stop("Nested categorical tests support A.method = \"delta\" only. ",
           "See `?semTests-support`.", call. = FALSE)
    }
    check_categorical_nested_pair(m0, m1)
    df <- as.integer(
      lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")
    )
    A <- get_a_matrix(m1, m0)
    if (nrow(A) != df) {
      stop("Nested categorical restriction rank (", nrow(A),
           ") does not match the df difference (", df, ").",
           call. = FALSE)
    }
    return(invisible(TRUE))
  }
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
