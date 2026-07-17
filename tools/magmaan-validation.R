#!/usr/bin/env Rscript

# Opt-in, non-CRAN parity checks against magmaan's independent C++23
# implementation. Run from the semTests package root:
#
#   Rscript tools/magmaan-validation.R
#
# This deliberately lives under tools/ instead of tests/testthat/: magmaan is
# not a CRAN package and should not become a package check dependency. pkgload
# is used so the checks exercise the current semTests working tree.

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to validate the current semTests source tree.",
       call. = FALSE)
}
if (!requireNamespace("magmaan", quietly = TRUE)) {
  stop(
    "Install magmaan's r-package/ directory before running this optional check.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
  library(lavaan)
})

spectrum_tolerance <- as.numeric(Sys.getenv(
  "SPECTRUM_TOLERANCE", "1e-4"
))
pvalue_tolerance <- as.numeric(Sys.getenv(
  "PVALUE_TOLERANCE", "1e-3"
))
if (!is.finite(spectrum_tolerance) || spectrum_tolerance <= 0 ||
    !is.finite(pvalue_tolerance) || pvalue_tolerance <= 0) {
  stop("Validation tolerances must be finite and positive.", call. = FALSE)
}

semtests_internal <- function(name) {
  get(name, envir = asNamespace("semTests"), inherits = FALSE)
}

check_close <- function(label, semtests, magmaan, tolerance) {
  semtests <- as.numeric(semtests)
  magmaan <- as.numeric(magmaan)
  if (length(semtests) != length(magmaan) ||
      any(!is.finite(semtests)) || any(!is.finite(magmaan))) {
    stop(label, ": non-finite values or unequal lengths.", call. = FALSE)
  }
  difference <- if (length(semtests)) {
    max(abs(semtests - magmaan))
  } else {
    0
  }
  if (difference > tolerance) {
    stop(
      label, ": maximum absolute difference ",
      format(difference, scientific = TRUE),
      " exceeds tolerance ", format(tolerance, scientific = TRUE), ".",
      call. = FALSE
    )
  }
  cat(sprintf(
    "%-46s max |difference| = %.3e\n",
    label, difference
  ))
  invisible(difference)
}

check_converged <- function(label, ...) {
  fits <- list(...)
  converged <- vapply(fits, function(fit) isTRUE(fit$converged), logical(1))
  if (!all(converged)) {
    stop(label, ": a magmaan fit did not converge.", call. = FALSE)
  }
}

tests <- c("SB", "SS", "PEBA4")

fiml_data <- function() {
  data <- lavaan::HolzingerSwineford1939
  for (j in seq_len(9L)) {
    data[
      seq(3L + j, nrow(data), by = 17L + j),
      paste0("x", j)
    ] <- NA_real_
  }
  data
}

fiml_h1 <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

fiml_h0 <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + a*x5 + a*x6
  speed   =~ x7 + x8 + x9
"

validate_fiml <- function() {
  data <- fiml_data()
  lavaan_h1 <- lavaan::cfa(
    fiml_h1, data, missing = "fiml", estimator = "MLR",
    meanstructure = TRUE
  )
  lavaan_h0 <- lavaan::cfa(
    fiml_h0, data, missing = "fiml", estimator = "MLR",
    meanstructure = TRUE
  )
  magmaan_h1 <- magmaan::magmaan(fiml_h1, data, estimator = "FIML")
  magmaan_h0 <- magmaan::magmaan(fiml_h0, data, estimator = "FIML")
  check_converged("FIML", magmaan_h1, magmaan_h0)

  df <- as.integer(lavaan::fitmeasures(lavaan_h1, "df"))
  semtests_spectrum <- semtests_internal("fiml_lambdas")(
    lavaan_h1, df, fiml.convention = "observed"
  )$ug_biased
  magmaan_tests <- magmaan::fmg_tests(magmaan_h1, tests = tests)
  check_close(
    "FIML single-model spectrum",
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance
  )
  check_close(
    "FIML single-model p-values",
    semTests::pvalues(
      lavaan_h1, tests, fiml.convention = "observed"
    ),
    magmaan_tests$p_value,
    pvalue_tolerance
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  for (A.method in c("delta", "exact")) {
    semtests_nested <- semtests_internal("fiml_lambdas_nested")(
      lavaan_h0, lavaan_h1, nested_df,
      A.method = A.method,
      fiml.convention = "observed"
    )$ug_biased
    magmaan_nested <- magmaan::fmg_nested(
      magmaan_h1, magmaan_h0,
      tests = tests, A.method = A.method
    )
    check_close(
      paste("FIML nested spectrum,", A.method),
      sort(semtests_nested),
      sort(magmaan_nested$eigenvalues[[1L]]),
      spectrum_tolerance
    )
    check_close(
      paste("FIML nested p-values,", A.method),
      semTests::pvalues_nested(
        lavaan_h0, lavaan_h1, tests = tests,
        A.method = A.method,
        fiml.convention = "observed"
      ),
      magmaan_nested$p_value,
      pvalue_tolerance
    )
  }
}

categorical_h1 <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

categorical_h0 <- "
  visual  =~ x1 + a*x2 + a*x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

categorical_data <- function(pairwise) {
  data <- lavaan::HolzingerSwineford1939
  ordered_names <- paste0("x", seq_len(9L))
  for (variable in ordered_names) {
    breaks <- c(
      -Inf,
      stats::quantile(data[[variable]], c(.33, .67)),
      Inf
    )
    data[[variable]] <- ordered(cut(
      data[[variable]],
      breaks = unique(breaks),
      include.lowest = TRUE
    ))
  }
  if (pairwise) {
    set.seed(314L)
    for (variable in ordered_names) {
      data[sample(nrow(data), 35L), variable] <- NA
    }
  }
  data
}

validate_categorical <- function(pairwise) {
  label <- if (pairwise) "pairwise" else "complete"
  missing <- if (pairwise) "pairwise" else "listwise"
  ordered_names <- paste0("x", seq_len(9L))
  data <- categorical_data(pairwise)

  lavaan_h1 <- lavaan::cfa(
    categorical_h1, data,
    ordered = ordered_names, estimator = "WLSMV", missing = missing
  )
  lavaan_h0 <- lavaan::cfa(
    categorical_h0, data,
    ordered = ordered_names, estimator = "WLSMV", missing = missing
  )

  magmaan_spec_h1 <- magmaan::model_spec(
    categorical_h1, ordered = ordered_names
  )
  magmaan_spec_h0 <- magmaan::model_spec(
    categorical_h0, ordered = ordered_names
  )
  # lavaan 0.7-2's pairwise Gamma uses observation-overlap denominators.
  magmaan_stats <- magmaan::magmaan_core$data_ordinal_stats_from_df(
    data, magmaan_spec_h1,
    ordered = ordered_names,
    missing = missing,
    pd_gamma = "overlap",
    full_wls_weight = FALSE
  )
  magmaan_h1 <- magmaan::magmaan(
    magmaan_spec_h1, magmaan_stats, estimator = "DWLS"
  )
  magmaan_h0 <- magmaan::magmaan(
    magmaan_spec_h0, magmaan_stats, estimator = "DWLS"
  )
  check_converged(paste("Categorical", label), magmaan_h1, magmaan_h0)

  df <- as.integer(lavaan::fitmeasures(lavaan_h1, "df"))
  semtests_spectrum <- semtests_internal("lavaan_lambdas")(
    lavaan_h1, df
  )$ug_biased
  magmaan_tests <- magmaan::fmg_tests_ordinal(
    magmaan_h1, magmaan_stats,
    tests = tests, weight = "DWLS"
  )
  check_close(
    paste("Categorical single spectrum,", label),
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance
  )
  check_close(
    paste("Categorical single p-values,", label),
    semTests::pvalues(lavaan_h1, tests),
    magmaan_tests$p_value,
    pvalue_tolerance
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  semtests_nested <- semtests_internal("lambdas_nested")(
    lavaan_h0, lavaan_h1,
    method = "2000", unbiased = 1L, df = nested_df
  )$ug_biased
  magmaan_nested <- magmaan::fmg_nested_ordinal(
    magmaan_h1, magmaan_h0, magmaan_stats,
    tests = tests, weight = "DWLS", A.method = "delta"
  )
  check_close(
    paste("Categorical nested spectrum,", label),
    sort(semtests_nested),
    sort(magmaan_nested$eigenvalues[[1L]]),
    spectrum_tolerance
  )
  check_close(
    paste("Categorical nested p-values,", label),
    semTests::pvalues_nested(
      lavaan_h0, lavaan_h1,
      tests = tests, method = "2000", A.method = "delta"
    ),
    magmaan_nested$p_value,
    pvalue_tolerance
  )
}

cat(sprintf(
  "semTests %s | lavaan %s | magmaan %s\n",
  as.character(utils::packageVersion("semTests")),
  as.character(utils::packageVersion("lavaan")),
  as.character(utils::packageVersion("magmaan"))
))
cat(sprintf(
  "Tolerances: spectrum %.1e, p-value %.1e\n\n",
  spectrum_tolerance, pvalue_tolerance
))

validate_fiml()
validate_categorical(pairwise = FALSE)
validate_categorical(pairwise = TRUE)

cat("\nAll optional magmaan parity checks passed.\n")
