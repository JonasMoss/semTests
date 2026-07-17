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
required_magmaan_exports <- c(
  "df_to_data", "fmg_tests", "fmg_nested_mixed_ordinal"
)
missing_magmaan_exports <- required_magmaan_exports[
  !vapply(
    required_magmaan_exports,
    exists, logical(1), envir = asNamespace("magmaan"), inherits = FALSE
  )
]
fmg_formals <- if (!length(missing_magmaan_exports)) {
  names(formals(get("fmg_tests", envir = asNamespace("magmaan"))))
} else {
  character()
}
if (length(missing_magmaan_exports) || !"gamma" %in% fmg_formals) {
  detail <- if (length(missing_magmaan_exports)) {
    paste0(" Missing: ", paste(missing_magmaan_exports, collapse = ", "), ".")
  } else {
    " `fmg_tests()` has no `gamma` argument."
  }
  stop(
    "The installed magmaan predates the parity R API used by this script.",
    detail,
    " Install the current magmaan r-package/ source and retry.",
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
endpoint_pvalue_tolerance <- as.numeric(Sys.getenv(
  "ENDPOINT_PVALUE_TOLERANCE", "2e-3"
))
multigroup_pvalue_tolerance <- as.numeric(Sys.getenv(
  "MULTIGROUP_PVALUE_TOLERANCE", "5e-3"
))
transform_tolerance <- as.numeric(Sys.getenv(
  "TRANSFORM_TOLERANCE", "2e-6"
))
statistic_tolerance <- as.numeric(Sys.getenv(
  "STATISTIC_TOLERANCE", "1e-4"
))
endpoint_statistic_tolerance <- as.numeric(Sys.getenv(
  "ENDPOINT_STATISTIC_TOLERANCE", "0.5"
))
parity_output <- Sys.getenv("PARITY_OUTPUT", "")
if (!is.finite(spectrum_tolerance) || spectrum_tolerance <= 0 ||
    !is.finite(pvalue_tolerance) || pvalue_tolerance <= 0 ||
    !is.finite(endpoint_pvalue_tolerance) ||
      endpoint_pvalue_tolerance <= 0 ||
    !is.finite(multigroup_pvalue_tolerance) ||
      multigroup_pvalue_tolerance <= 0 ||
    !is.finite(transform_tolerance) || transform_tolerance <= 0 ||
    !is.finite(statistic_tolerance) || statistic_tolerance <= 0 ||
    !is.finite(endpoint_statistic_tolerance) ||
      endpoint_statistic_tolerance <= 0) {
  stop("Validation tolerances must be finite and positive.", call. = FALSE)
}

semtests_internal <- function(name) {
  get(name, envir = asNamespace("semTests"), inherits = FALSE)
}

parity_results <- list()

check_close <- function(label, semtests, magmaan, tolerance,
                        layer = "end-to-end") {
  value_names <- names(semtests)
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
  parity_results[[length(parity_results) + 1L]] <<- data.frame(
    layer = layer,
    check = label,
    values = length(semtests),
    max_abs_difference = difference,
    tolerance = tolerance,
    stringsAsFactors = FALSE
  )
  if (difference > tolerance) {
    worst <- which.max(abs(semtests - magmaan))
    where <- if (length(value_names) >= worst && nzchar(value_names[worst])) {
      paste0(" at `", value_names[worst], "`")
    } else {
      paste0(" at value ", worst)
    }
    stop(
      label, ": maximum absolute difference ",
      format(difference, scientific = TRUE), where,
      " (semTests = ", format(semtests[worst], digits = 12L),
      ", magmaan = ", format(magmaan[worst], digits = 12L), ")",
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

check_names <- function(label, semtests, magmaan) {
  semtests <- as.character(semtests)
  magmaan <- as.character(magmaan)
  if (!identical(semtests, magmaan)) {
    stop(
      label, ": labels differ.\nsemTests: ", paste(semtests, collapse = ", "),
      "\nmagmaan: ", paste(magmaan, collapse = ", "),
      call. = FALSE
    )
  }
  cat(sprintf("%-46s exact label match (%d)\n", label, length(semtests)))
  invisible(TRUE)
}

check_converged <- function(label, ...) {
  fits <- list(...)
  converged <- vapply(fits, function(fit) isTRUE(fit$converged), logical(1))
  if (!all(converged)) {
    stop(label, ": a magmaan fit did not converge.", call. = FALSE)
  }
}

tests_for_df <- function(df, base = NULL, unbiased = FALSE) {
  if (length(df) != 1L || !is.finite(df) || df < 1L) {
    stop("Validation requires a positive df.", call. = FALSE)
  }
  blocks <- min(4L, as.integer(df))
  tests <- c(
    "STD", "SB", "SS", "SF", "ALL", "PALL",
    paste0("EBA", blocks), paste0("PEBA", blocks),
    "POLS0.5", "POLS2", "POLS6"
  )
  if (unbiased) {
    tests <- ifelse(tests == "STD", tests, paste0(tests, "_UG"))
  }
  if (!is.null(base)) {
    tests <- paste0(tests, "_", toupper(base))
  }
  tests
}

validate_transforms <- function() {
  cases <- list(
    five_df = list(
      statistic = 8.75,
      spectrum = c(4.6, 3.2, 2.1, 1.4, 0.9)
    ),
    uneven_blocks = list(
      statistic = 17.5,
      spectrum = c(2.4, 1.9, 1.5, 1.1, 0.75, 0.4, 0.18)
    ),
    concentrated = list(
      statistic = 6.25,
      spectrum = c(3.8, 0.9, 0.2)
    )
  )
  methods <- list(
    standard = list(magmaan = "standard", semtests = "std", param = 4),
    sb = list(magmaan = "sb", semtests = "sb", param = 4),
    ss = list(magmaan = "ss", semtests = "ss", param = 4),
    sf = list(magmaan = "scaled_f", semtests = "sf", param = 4),
    all = list(magmaan = "all", semtests = "all", param = 4),
    pall = list(magmaan = "penalized_all", semtests = "pall", param = 4),
    eba = list(magmaan = "eba", semtests = "eba", param = 2),
    peba = list(magmaan = "peba", semtests = "peba", param = 2),
    pols_half = list(magmaan = "pols", semtests = "pols", param = 0.5),
    pols_two = list(magmaan = "pols", semtests = "pols", param = 2),
    pols_six = list(magmaan = "pols", semtests = "pols", param = 6)
  )

  semtests_pvalue <- function(method, statistic, spectrum, param) {
    switch(
      method,
      std = semtests_internal("trad_pvalue")(
        length(spectrum), statistic, spectrum, "std"
      ),
      sb = semtests_internal("trad_pvalue")(
        length(spectrum), statistic, spectrum, "sb"
      ),
      ss = semtests_internal("trad_pvalue")(
        length(spectrum), statistic, spectrum, "ss"
      ),
      sf = semtests_internal("trad_pvalue")(
        length(spectrum), statistic, spectrum, "sf"
      ),
      all = semtests_internal("trad_pvalue")(
        length(spectrum), statistic, spectrum, "all"
      ),
      pall = semtests_internal("trad_pvalue")(
        length(spectrum), statistic, spectrum, "pall"
      ),
      eba = semtests_internal("eba_pvalue")(statistic, spectrum, param),
      peba = semtests_internal("peba_pvalue")(statistic, spectrum, param),
      pols = semtests_internal("pols_pvalue")(statistic, spectrum, param)
    )
  }

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    semtests <- magmaan <- numeric(length(methods))
    for (i in seq_along(methods)) {
      method <- methods[[i]]
      semtests[i] <- semtests_pvalue(
        method$semtests, case$statistic, case$spectrum, method$param
      )
      magmaan[i] <- magmaan::magmaan_core$robust_fmg_test(
        case$statistic, length(case$spectrum), case$spectrum,
        method = method$magmaan, param = method$param
      )$p_value
    }
    check_close(
      paste("Fixed-spectrum transforms,", case_name),
      semtests, magmaan, transform_tolerance, layer = "transform"
    )
  }
}

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
  visual  =~ x1 + b*x2 + b*x3
  textual =~ x4 + a*x5 + a*x6
  speed   =~ x7 + x8 + x9
"

validate_classical_ml <- function() {
  data <- lavaan::HolzingerSwineford1939
  lavaan_h1 <- lavaan::cfa(fiml_h1, data, estimator = "MLM")
  lavaan_h0 <- lavaan::cfa(fiml_h0, data, estimator = "MLM")
  magmaan_h1 <- magmaan::magmaan(fiml_h1, data, estimator = "ML")
  magmaan_h0 <- magmaan::magmaan(fiml_h0, data, estimator = "ML")
  check_converged("Classical ML", magmaan_h1, magmaan_h0)

  df <- as.integer(lavaan::fitmeasures(lavaan_h1, "df"))
  single_tests <- unique(c(
    tests_for_df(df, "ml"),
    tests_for_df(df, "rls"),
    tests_for_df(df, "ml", unbiased = TRUE),
    tests_for_df(df, "rls", unbiased = TRUE)
  ))
  semtests_single <- semTests::pvalues(lavaan_h1, single_tests)
  magmaan_single <- magmaan::fmg_tests(magmaan_h1, tests = single_tests)
  check_names(
    "Classical single-model labels",
    names(semtests_single), magmaan_single$label
  )

  semtests_spectra <- lapply(
    semtests_internal("ugamma")(lavaan_h1, unbiased = 3L),
    semtests_internal("ugamma_eigenvalues"),
    df = df
  )
  biased_row <- which(!magmaan_single$ug)[1L]
  unbiased_row <- which(magmaan_single$ug)[1L]
  check_close(
    "Classical single spectrum, biased",
    sort(semtests_spectra$ug_biased),
    sort(magmaan_single$eigenvalues[[biased_row]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    "Classical single spectrum, unbiased",
    sort(semtests_spectra$ug_unbiased),
    sort(magmaan_single$eigenvalues[[unbiased_row]]),
    spectrum_tolerance, layer = "spectrum"
  )

  semtests_statistics <- semtests_internal("make_chisqs")(
    c("ml", "rls"), lavaan_h1
  )
  magmaan_statistics <- vapply(
    names(semtests_statistics),
    function(base) {
      magmaan_single$base_statistic[match(base, magmaan_single$base)]
    },
    numeric(1)
  )
  check_close(
    "Classical single base statistics",
    semtests_statistics, magmaan_statistics,
    statistic_tolerance, layer = "statistic"
  )
  check_close(
    "Classical single-model p-values",
    semtests_single, magmaan_single$p_value,
    pvalue_tolerance, layer = "end-to-end"
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  nested_tests <- c(
    tests_for_df(nested_df, "ml"),
    tests_for_df(nested_df, "rls")
  )
  semtests_nested <- semTests::pvalues_nested(
    lavaan_h0, lavaan_h1,
    method = "2000", tests = nested_tests
  )
  # semTests' complete-data method-2000 construction is the local
  # moment-tangent map. magmaan must use `delta` explicitly; its `exact`
  # default is a different restriction construction.
  magmaan_nested <- magmaan::fmg_nested(
    magmaan_h1, magmaan_h0, data = data,
    tests = nested_tests, A.method = "delta"
  )
  check_names(
    "Classical nested labels",
    names(semtests_nested), magmaan_nested$label
  )
  semtests_nested_spectrum <- semtests_internal("lambdas_nested")(
    lavaan_h0, lavaan_h1,
    method = "2000", unbiased = 1L, df = nested_df
  )$ug_biased
  check_close(
    "Classical nested spectrum, method 2000 delta",
    sort(semtests_nested_spectrum),
    sort(magmaan_nested$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  semtests_nested_statistics <- semtests_internal("make_chisqs")(
    c("ml", "rls"), lavaan_h0, lavaan_h1
  )
  magmaan_nested_statistics <- vapply(
    names(semtests_nested_statistics),
    function(base) {
      magmaan_nested$base_statistic[match(base, magmaan_nested$base)]
    },
    numeric(1)
  )
  check_close(
    "Classical nested base statistics",
    semtests_nested_statistics, magmaan_nested_statistics,
    statistic_tolerance, layer = "statistic"
  )
  check_close(
    "Classical nested p-values, method 2000 delta",
    semtests_nested, magmaan_nested$p_value,
    pvalue_tolerance, layer = "end-to-end"
  )
}

validate_continuous_ls <- function(estimator) {
  data <- lavaan::HolzingerSwineford1939
  lavaan_args <- list(
    model = fiml_h1, data = data, estimator = estimator
  )
  if (identical(estimator, "ULS")) {
    lavaan_args$test <- "satorra.bentler"
  }
  lavaan_fit <- do.call(lavaan::cfa, lavaan_args)
  magmaan_spec <- magmaan::model_spec(fiml_h1)
  # lavaan's continuous LS sample covariance uses the N-1 convention. Magmaan
  # defaults to N for data-frame fits, so request the matching public data
  # convention explicitly instead of treating the scaling difference as an
  # implementation discrepancy.
  magmaan_data <- magmaan::df_to_data(
    data, magmaan_spec, scaling = "n-1"
  )
  magmaan_fit <- magmaan::magmaan(
    magmaan_spec, magmaan_data, estimator = estimator
  )
  check_converged(paste("Continuous", estimator), magmaan_fit)

  df <- as.integer(lavaan::fitmeasures(lavaan_fit, "df"))
  tests <- tests_for_df(df)
  semtests_spectrum <- semtests_internal("lavaan_lambdas")(
    lavaan_fit, df
  )$ug_biased
  # lavaan's ULS family defaults to robust.sem.nt (normal-theory Gamma);
  # GLS uses empirical Gamma. Exercise both public magmaan routes.
  gamma <- if (identical(estimator, "ULS")) "normal" else "empirical"
  magmaan_tests <- magmaan::fmg_tests(
    magmaan_fit, tests = tests, gamma = gamma
  )
  check_close(
    paste("Continuous single spectrum,", estimator),
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    paste("Continuous single p-values,", estimator),
    semTests::pvalues(lavaan_fit, tests),
    magmaan_tests$p_value,
    endpoint_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    paste("Continuous single base statistic,", estimator),
    semtests_internal("make_chisqs")("ml", lavaan_fit),
    magmaan_tests$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )
}

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
  single_tests <- tests_for_df(df)
  semtests_spectrum <- semtests_internal("fiml_lambdas")(
    lavaan_h1, df, fiml.convention = "observed"
  )$ug_biased
  magmaan_tests <- magmaan::fmg_tests(magmaan_h1, tests = single_tests)
  check_close(
    "FIML single-model spectrum",
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  semtests_pvalues <- semTests::pvalues(
    lavaan_h1, single_tests, fiml.convention = "observed"
  )
  check_names(
    "FIML single-model labels",
    names(semtests_pvalues), magmaan_tests$label
  )
  check_close(
    "FIML single-model p-values",
    semtests_pvalues,
    magmaan_tests$p_value,
    pvalue_tolerance, layer = "end-to-end"
  )
  check_close(
    "FIML single-model base statistic",
    semtests_internal("make_chisqs")("ml", lavaan_h1),
    magmaan_tests$base_statistic[1L],
    statistic_tolerance, layer = "statistic"
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  nested_tests <- tests_for_df(nested_df)
  for (A.method in c("delta", "exact")) {
    semtests_nested <- semtests_internal("fiml_lambdas_nested")(
      lavaan_h0, lavaan_h1, nested_df,
      A.method = A.method,
      fiml.convention = "observed"
    )$ug_biased
    magmaan_nested <- magmaan::fmg_nested(
      magmaan_h1, magmaan_h0,
      tests = nested_tests, A.method = A.method
    )
    check_close(
      paste("FIML nested spectrum,", A.method),
      sort(semtests_nested),
      sort(magmaan_nested$eigenvalues[[1L]]),
      spectrum_tolerance, layer = "spectrum"
    )
    semtests_pvalues <- semTests::pvalues_nested(
      lavaan_h0, lavaan_h1, tests = nested_tests,
      A.method = A.method,
      fiml.convention = "observed"
    )
    check_names(
      paste("FIML nested labels,", A.method),
      names(semtests_pvalues), magmaan_nested$label
    )
    check_close(
      paste("FIML nested p-values,", A.method),
      semtests_pvalues,
      magmaan_nested$p_value,
      pvalue_tolerance, layer = "end-to-end"
    )
    check_close(
      paste("FIML nested base statistic,", A.method),
      semtests_internal("make_chisqs")("ml", lavaan_h0, lavaan_h1),
      magmaan_nested$base_statistic[1L],
      statistic_tolerance, layer = "statistic"
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

validate_categorical <- function(pairwise,
                                 lavaan_estimator = "WLSMV",
                                 magmaan_estimator = "DWLS") {
  label <- paste(
    magmaan_estimator,
    if (pairwise) "pairwise" else "complete"
  )
  missing <- if (pairwise) "pairwise" else "listwise"
  ordered_names <- paste0("x", seq_len(9L))
  data <- categorical_data(pairwise)

  lavaan_h1 <- lavaan::cfa(
    categorical_h1, data,
    ordered = ordered_names,
    estimator = lavaan_estimator, missing = missing
  )
  lavaan_h0 <- lavaan::cfa(
    categorical_h0, data,
    ordered = ordered_names,
    estimator = lavaan_estimator, missing = missing
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
    full_wls_weight = identical(magmaan_estimator, "WLS")
  )
  magmaan_h1 <- magmaan::magmaan(
    magmaan_spec_h1, magmaan_stats, estimator = magmaan_estimator
  )
  magmaan_h0 <- magmaan::magmaan(
    magmaan_spec_h0, magmaan_stats, estimator = magmaan_estimator
  )
  check_converged(paste("Categorical", label), magmaan_h1, magmaan_h0)

  df <- as.integer(lavaan::fitmeasures(lavaan_h1, "df"))
  single_tests <- tests_for_df(df)
  semtests_spectrum <- semtests_internal("lavaan_lambdas")(
    lavaan_h1, df
  )$ug_biased
  magmaan_tests <- magmaan::fmg_tests_ordinal(
    magmaan_h1, magmaan_stats,
    tests = single_tests, weight = magmaan_estimator
  )
  check_close(
    paste("Categorical single spectrum,", label),
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    paste("Categorical single p-values,", label),
    semTests::pvalues(lavaan_h1, single_tests),
    magmaan_tests$p_value,
    endpoint_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    paste("Categorical single base statistic,", label),
    semtests_internal("make_chisqs")("ml", lavaan_h1),
    magmaan_tests$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  nested_tests <- tests_for_df(nested_df)
  semtests_nested <- semtests_internal("lambdas_nested")(
    lavaan_h0, lavaan_h1,
    method = "2000", unbiased = 1L, df = nested_df
  )$ug_biased
  magmaan_nested <- magmaan::fmg_nested_ordinal(
    magmaan_h1, magmaan_h0, magmaan_stats,
    tests = nested_tests, weight = magmaan_estimator, A.method = "delta"
  )
  check_close(
    paste("Categorical nested spectrum,", label),
    sort(semtests_nested),
    sort(magmaan_nested$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    paste("Categorical nested p-values,", label),
    semTests::pvalues_nested(
      lavaan_h0, lavaan_h1,
      tests = nested_tests, method = "2000", A.method = "delta"
    ),
    magmaan_nested$p_value,
    endpoint_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    paste("Categorical nested base statistic,", label),
    semtests_internal("make_chisqs")("ml", lavaan_h0, lavaan_h1),
    magmaan_nested$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )
}

validate_mixed_categorical <- function() {
  data <- lavaan::HolzingerSwineford1939
  ordered_names <- paste0("x", seq_len(3L))
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

  lavaan_h1 <- lavaan::cfa(
    categorical_h1, data,
    ordered = ordered_names, estimator = "WLSMV",
    meanstructure = TRUE
  )
  lavaan_h0 <- lavaan::cfa(
    categorical_h0, data,
    ordered = ordered_names, estimator = "WLSMV",
    meanstructure = TRUE
  )
  magmaan_spec_h1 <- magmaan::model_spec(
    categorical_h1,
    ordered = ordered_names, meanstructure = TRUE
  )
  magmaan_spec_h0 <- magmaan::model_spec(
    categorical_h0,
    ordered = ordered_names, meanstructure = TRUE
  )
  magmaan_stats <- magmaan::magmaan_core$data_mixed_ordinal_stats_from_df(
    data, magmaan_spec_h1,
    ordered = ordered_names, missing = "listwise",
    full_wls_weight = FALSE
  )
  magmaan_h1 <- magmaan::magmaan(
    magmaan_spec_h1, magmaan_stats, estimator = "DWLS"
  )
  magmaan_h0 <- magmaan::magmaan(
    magmaan_spec_h0, magmaan_stats, estimator = "DWLS"
  )
  check_converged("Mixed categorical DWLS", magmaan_h1, magmaan_h0)

  df <- as.integer(lavaan::fitmeasures(lavaan_h1, "df"))
  tests <- tests_for_df(df)
  semtests_spectrum <- semtests_internal("lavaan_lambdas")(
    lavaan_h1, df
  )$ug_biased
  magmaan_tests <- magmaan::fmg_tests_mixed_ordinal(
    magmaan_h1, magmaan_stats,
    tests = tests, weight = "DWLS"
  )
  check_close(
    "Mixed categorical single spectrum, DWLS",
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    "Mixed categorical single p-values, DWLS",
    semTests::pvalues(lavaan_h1, tests),
    magmaan_tests$p_value,
    endpoint_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    "Mixed categorical base statistic, DWLS",
    semtests_internal("make_chisqs")("ml", lavaan_h1),
    magmaan_tests$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  nested_tests <- tests_for_df(nested_df)
  semtests_nested_spectrum <- semtests_internal("lambdas_nested")(
    lavaan_h0, lavaan_h1,
    method = "2000", unbiased = 1L, df = nested_df
  )$ug_biased
  magmaan_nested <- magmaan::fmg_nested_mixed_ordinal(
    magmaan_h1, magmaan_h0, magmaan_stats,
    tests = nested_tests, weight = "DWLS", A.method = "delta"
  )
  check_close(
    "Mixed categorical nested spectrum, DWLS",
    sort(semtests_nested_spectrum),
    sort(magmaan_nested$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    "Mixed categorical nested p-values, DWLS",
    semTests::pvalues_nested(
      lavaan_h0, lavaan_h1,
      tests = nested_tests, method = "2000", A.method = "delta"
    ),
    magmaan_nested$p_value,
    endpoint_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    "Mixed categorical nested statistic, DWLS",
    semtests_internal("make_chisqs")("ml", lavaan_h0, lavaan_h1),
    magmaan_nested$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )
}

validate_multigroup_categorical <- function() {
  data <- categorical_data(pairwise = FALSE)
  ordered_names <- paste0("x", seq_len(9L))
  group_labels <- unique(as.character(data$school))

  lavaan_h1 <- lavaan::cfa(
    categorical_h1, data,
    ordered = ordered_names, group = "school",
    estimator = "WLSMV"
  )
  lavaan_h0 <- lavaan::cfa(
    categorical_h1, data,
    ordered = ordered_names, group = "school",
    group.equal = "loadings", estimator = "WLSMV"
  )
  magmaan_spec_h1 <- magmaan::model_spec(
    categorical_h1,
    ordered = ordered_names,
    group = "school", group_labels = group_labels
  )
  magmaan_spec_h0 <- magmaan::model_spec(
    categorical_h1,
    ordered = ordered_names,
    group = "school", group_labels = group_labels,
    group_equal = "loadings"
  )
  magmaan_stats <- magmaan::magmaan_core$data_ordinal_stats_from_df(
    data, magmaan_spec_h1,
    ordered = ordered_names, group = "school",
    missing = "listwise", pd_gamma = "overlap",
    full_wls_weight = FALSE
  )
  magmaan_h1 <- magmaan::magmaan(
    magmaan_spec_h1, magmaan_stats, estimator = "DWLS"
  )
  magmaan_h0 <- magmaan::magmaan(
    magmaan_spec_h0, magmaan_stats, estimator = "DWLS"
  )
  check_converged("Multigroup categorical DWLS", magmaan_h1, magmaan_h0)

  df <- as.integer(lavaan::fitmeasures(lavaan_h1, "df"))
  single_tests <- tests_for_df(df)
  semtests_spectrum <- semtests_internal("lavaan_lambdas")(
    lavaan_h1, df
  )$ug_biased
  magmaan_tests <- magmaan::fmg_tests_ordinal(
    magmaan_h1, magmaan_stats,
    tests = single_tests, weight = "DWLS"
  )
  check_close(
    "Multigroup categorical single spectrum, DWLS",
    sort(semtests_spectrum),
    sort(magmaan_tests$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    "Multigroup categorical single p-values, DWLS",
    semTests::pvalues(lavaan_h1, single_tests),
    magmaan_tests$p_value,
    multigroup_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    "Multigroup categorical single statistic, DWLS",
    semtests_internal("make_chisqs")("ml", lavaan_h1),
    magmaan_tests$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )

  nested_df <- as.integer(
    lavaan::fitmeasures(lavaan_h0, "df") -
      lavaan::fitmeasures(lavaan_h1, "df")
  )
  nested_tests <- tests_for_df(nested_df)
  semtests_nested_spectrum <- semtests_internal("lambdas_nested")(
    lavaan_h0, lavaan_h1,
    method = "2000", unbiased = 1L, df = nested_df
  )$ug_biased
  magmaan_nested <- magmaan::fmg_nested_ordinal(
    magmaan_h1, magmaan_h0, magmaan_stats,
    tests = nested_tests, weight = "DWLS", A.method = "delta"
  )
  check_close(
    "Multigroup categorical nested spectrum, DWLS",
    sort(semtests_nested_spectrum),
    sort(magmaan_nested$eigenvalues[[1L]]),
    spectrum_tolerance, layer = "spectrum"
  )
  check_close(
    "Multigroup categorical nested p-values, DWLS",
    semTests::pvalues_nested(
      lavaan_h0, lavaan_h1,
      tests = nested_tests, method = "2000", A.method = "delta"
    ),
    magmaan_nested$p_value,
    multigroup_pvalue_tolerance, layer = "optimizer endpoint"
  )
  check_close(
    "Multigroup categorical nested statistic, DWLS",
    semtests_internal("make_chisqs")("ml", lavaan_h0, lavaan_h1),
    magmaan_nested$base_statistic[1L],
    endpoint_statistic_tolerance, layer = "optimizer endpoint"
  )
}

package_fingerprint <- function(package) {
  package_dir <- system.file(package = package)
  candidates <- c(
    file.path(
      package_dir, "libs",
      paste0(package, .Platform$dynlib.ext)
    ),
    file.path(package_dir, "R", paste0(package, ".rdb"))
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return("unavailable")
  paste(unname(tools::md5sum(candidates)), collapse = ":")
}

semtests_version <- as.character(utils::packageVersion("semTests"))
lavaan_version <- as.character(utils::packageVersion("lavaan"))
magmaan_version <- as.character(utils::packageVersion("magmaan"))
magmaan_fingerprint <- package_fingerprint("magmaan")

cat(sprintf(
  "semTests %s | lavaan %s | magmaan %s\n",
  semtests_version, lavaan_version, magmaan_version
))
cat("magmaan installed-code fingerprint:", magmaan_fingerprint, "\n")
cat(sprintf(
  paste0(
    "Tolerances: transform %.1e, statistic %.1e, ",
    "spectrum %.1e, p-value %.1e, optimizer p-value %.1e, ",
    "multigroup p-value %.1e, optimizer statistic %.1e\n\n"
  ),
  transform_tolerance, statistic_tolerance,
  spectrum_tolerance, pvalue_tolerance,
  endpoint_pvalue_tolerance, multigroup_pvalue_tolerance,
  endpoint_statistic_tolerance
))

validate_transforms()
validate_classical_ml()
validate_continuous_ls("GLS")
validate_continuous_ls("ULS")
validate_fiml()
validate_categorical(
  pairwise = FALSE,
  lavaan_estimator = "WLSMV", magmaan_estimator = "DWLS"
)
validate_categorical(
  pairwise = TRUE,
  lavaan_estimator = "WLSMV", magmaan_estimator = "DWLS"
)
validate_categorical(
  pairwise = FALSE,
  lavaan_estimator = "ULSMV", magmaan_estimator = "ULS"
)
validate_categorical(
  pairwise = TRUE,
  lavaan_estimator = "ULSMV", magmaan_estimator = "ULS"
)
validate_categorical(
  pairwise = FALSE,
  lavaan_estimator = "WLS", magmaan_estimator = "WLS"
)
validate_categorical(
  pairwise = TRUE,
  lavaan_estimator = "WLS", magmaan_estimator = "WLS"
)
validate_mixed_categorical()
validate_multigroup_categorical()

parity_results <- do.call(rbind, parity_results)
parity_results$semTests_version <- semtests_version
parity_results$lavaan_version <- lavaan_version
parity_results$magmaan_version <- magmaan_version
parity_results$magmaan_fingerprint <- magmaan_fingerprint

if (nzchar(parity_output)) {
  dir.create(dirname(parity_output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(parity_results, parity_output, row.names = FALSE)
  cat("Machine-readable parity report:", parity_output, "\n")
}

cat(sprintf(
  "\nAll %d optional magmaan parity checks passed (%d compared values).\n",
  nrow(parity_results), sum(parity_results$values)
))
