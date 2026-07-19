# Does the fitted model contain observed exogenous predictors?
has_observed_exogenous <- function(fit) {
  any(fit@Model@nexo > 0L)
}

# Is inference conditional on observed exogenous predictors?
uses_fixed_or_conditional_x <- function(fit) {
  has_observed_exogenous(fit) &&
    (isTRUE(fit@Options$fixed.x) ||
      isTRUE(fit@Options$conditional.x) ||
      isTRUE(fit@Model@conditional.x))
}

#' Is this the classical normal-theory, complete-data, random-x ML case?
#'
#' The Du-Bentler unbiased gamma and the RLS (`browne.residual.nt.model`)
#' statistic are defined here. Other families use the biased gamma and the
#' estimator's own uncorrected statistic.
#' @keywords internal
is_classic_nt <- function(fit) {
  fit@Options$estimator == "ML" &&
    !isTRUE(fit@Model@categorical) &&
    identical(fit@Options$missing, "listwise") &&
    !uses_fixed_or_conditional_x(fit)
}

#' @keywords internal
make_chisqs <- function(chisq, m0, m1) {
  # `lavTest()` returns one of two shapes depending on the requested test.
  # Asking for "standard" gives a single flat test object with a top-level
  # `$stat`; asking for a non-standard test (e.g. "browne.residual.nt.model")
  # gives a named list keyed by "standard" and the requested test, each element
  # carrying its own `$stat`. Read the statistic from whichever shape we get.
  get_stat <- function(object, test) {
    res <- lavaan::lavTest(object, test = test)
    if (!is.null(res[[test]])) res[[test]][["stat"]] else res[["stat"]]
  }
  # "ml" is the standard (uncorrected) test statistic of the fit -- the LRT for
  # ML, and the corresponding uncorrected quadratic form for ULS/GLS/DWLS/FIML.
  # It is the right input to the eigenvalue correction for any estimator; the
  # native scaled/shifted statistic is already corrected and must NOT be used.
  standard <- function(object) get_stat(object, "standard")
  rls <- function(object) get_stat(object, "browne.residual.nt.model")
  wrap <- function(f, object) if (missing(object)) 0 else f(object)

  classic_nt <- is_classic_nt(m0)
  # The default statistic is fit-appropriate: RLS for the classical case (its
  # historical default), the standard statistic otherwise.
  chisq[chisq == "auto"] <- if (classic_nt) "rls" else "ml"
  chisq <- unique(chisq)

  if ("rls" %in% chisq && !classic_nt) {
    stop("The RLS statistic ('browne.residual.nt.model') is only available for ",
      "continuous, complete-data ML with random observed exogenous ",
      "covariates (`fixed.x = FALSE`, `conditional.x = FALSE`). lavaan ",
      "evaluates it as ADF outside that setting. ",
      "Use the standard statistic (the `_ML` suffix) or omit the suffix. ",
      "See `?semTests-support`.",
      call. = FALSE
    )
  }

  chisqs <- c()
  if ("ml" %in% chisq) chisqs["ml"] <- standard(m0) - wrap(standard, m1)
  if ("rls" %in% chisq) chisqs["rls"] <- rls(m0) - wrap(rls, m1)
  chisqs
}

#' Extract model degrees of freedom through lavaan's public API.
#' @keywords internal
fit_df <- function(fit) {
  unname(lavaan::fitmeasures(fit, "df"))
}

#' Parse and validate public test specifications.
#' @keywords internal
parse_tests <- function(tests) {
  if (is.null(tests)) {
    semtests_abort(
      "`tests` must request at least one p-value; see `?pvalues` for valid names.",
      "semTests_error_invalid_tests"
    )
  }
  validate_tests(tests)
  lapply(tests, split_input)
}

#' Validate public test specifications.
#'
#' @param tests Character vector supplied to the public `tests` argument.
#' @return `tests`, invisibly.
#' @keywords internal
validate_tests <- function(tests) {
  if (!is.character(tests) || !length(tests) || anyNA(tests) ||
    any(!nzchar(tests))) {
    semtests_abort(
      "`tests` must be a non-empty character vector with no missing or empty values.",
      "semTests_error_invalid_tests"
    )
  }
  invisible(tests)
}

#' Split string into options.
#' @param string Input string
#' @keywords internal

split_input <- function(string) {
  validate_tests(string)
  if (length(string) != 1L) {
    semtests_abort(
      "`split_input()` requires exactly one test specification.",
      "semTests_error_invalid_tests"
    )
  }
  string <- tolower(string)
  splitted <- strsplit(string, "_", fixed = TRUE)[[1]]
  valid_suffix <- length(splitted) == 1L ||
    (length(splitted) == 2L &&
      splitted[2L] %in% c("ug", "ml", "rls")) ||
    (length(splitted) == 3L &&
      splitted[2L] == "ug" &&
      splitted[3L] %in% c("ml", "rls"))
  if (!valid_suffix) {
    semtests_abort(
      paste0(
        "Invalid test specification `", string, "`. Use ",
        "`TEST`, `TEST_UG`, `TEST_ML`, `TEST_RLS`, `TEST_UG_ML`, or ",
        "`TEST_UG_RLS`; see `?pvalues`."
      ),
      "semTests_error_invalid_tests"
    )
  }
  trad <- peba <- eba <- pols <- NULL
  unbiased <- 1
  chisq <- "auto" # resolved per fit in make_chisqs (rls if classical NT, else ml)
  type <- splitted[1L]

  if (length(splitted) == 3) {
    unbiased <- if (splitted[2] == "ug") 2 else 1
    chisq <- splitted[3]
  } else if (length(splitted) == 2) {
    if (splitted[2] == "rls" || splitted[2] == "ml") {
      chisq <- splitted[2]
    } else if (splitted[2] == "ug") {
      unbiased <- 2
    }
  }

  numeric_parameter <- function(prefix, integer = FALSE) {
    value_string <- substring(type, nchar(prefix) + 1L)
    value <- if (value_string == "") 2 else suppressWarnings(as.numeric(value_string))
    valid <- length(value) == 1L && is.finite(value) && value > 0
    if (integer) {
      valid <- valid && value == floor(value) && value <= .Machine$integer.max
    }
    if (!valid) {
      requirement <- if (integer) "a positive integer" else "a finite positive number"
      semtests_abort(
        paste0(
          "Invalid test specification `", string, "`: `", prefix,
          "` requires ", requirement, "."
        ),
        "semTests_error_invalid_tests"
      )
    }
    if (integer) as.integer(value) else value
  }

  if (startsWith(type, "peba")) {
    peba <- numeric_parameter("peba", integer = TRUE)
  } else if (startsWith(type, "eba")) {
    eba <- numeric_parameter("eba", integer = TRUE)
  } else if (startsWith(type, "pols")) {
    pols <- numeric_parameter("pols")
  } else if (type %in% c("std", "sb", "sf", "ss", "all", "pall")) {
    trad <- type
  } else {
    semtests_abort(
      paste0(
        "Unknown test family in `", string,
        "`. See `?pvalues` for the supported test names."
      ),
      "semTests_error_invalid_tests"
    )
  }

  list(
    trad = trad,
    eba = eba,
    peba = peba,
    pols = pols,
    unbiased = unbiased,
    chisq = chisq
  )
}
