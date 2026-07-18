categorical_test_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

categorical_test_model0 <- "
  visual  =~ x1 + a*x2 + a*x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

categorical_test_data <- function(binary = FALSE, mixed = FALSE,
                                  pairwise = FALSE, seed = 314L) {
  data <- lavaan::HolzingerSwineford1939
  ordered <- if (mixed) paste0("x", 1:6) else paste0("x", 1:9)
  for (variable in ordered) {
    breaks <- if (binary) {
      c(-Inf, stats::median(data[[variable]]), Inf)
    } else {
      c(-Inf, stats::quantile(data[[variable]], c(.33, .67)), Inf)
    }
    data[[variable]] <- ordered(cut(
      data[[variable]], breaks = unique(breaks),
      include.lowest = TRUE
    ))
  }
  if (pairwise) {
    set.seed(seed)
    for (variable in paste0("x", 1:9)) {
      data[sample(nrow(data), 35L), variable] <- NA
    }
  }
  list(data = data, ordered = ordered)
}

fit_categorical_pair <- function(estimator = "WLSMV", pairwise = FALSE,
                                 mixed = FALSE, group = NULL,
                                 parameterization = "delta") {
  fixture <- categorical_test_data(pairwise = pairwise, mixed = mixed)
  arguments <- list(
    data = fixture$data,
    ordered = fixture$ordered,
    estimator = estimator,
    missing = if (pairwise) "pairwise" else "listwise",
    group = group,
    parameterization = parameterization
  )
  list(
    m0 = do.call(lavaan::cfa, c(list(categorical_test_model0), arguments)),
    m1 = do.call(lavaan::cfa, c(list(categorical_test_model), arguments))
  )
}

categorical_valid_p <- function(x) {
  all(is.finite(x) & x >= 0 & x <= 1)
}

categorical_fixture <- categorical_test_data()
categorical_base_estimators <- c(
  "WLSMV", "DWLS", "WLSM", "WLSMVS", "ULSMV", "ULS"
)
categorical_base_fits <- setNames(
  lapply(categorical_base_estimators, function(estimator) {
    lavaan::cfa(
      categorical_test_model, categorical_fixture$data,
      ordered = categorical_fixture$ordered, estimator = estimator
    )
  }),
  categorical_base_estimators
)
categorical_pairs <- list(
  default = fit_categorical_pair(),
  pairwise = fit_categorical_pair(pairwise = TRUE),
  mixed = fit_categorical_pair(mixed = TRUE),
  theta = fit_categorical_pair(parameterization = "theta")
)
categorical_multigroup <- lavaan::cfa(
  categorical_test_model, categorical_fixture$data,
  ordered = categorical_fixture$ordered, group = "school"
)

test_that("biased single-model spectra come directly from lavaan UGamma", {
  for (estimator in names(categorical_base_fits)) {
    fit <- categorical_base_fits[[estimator]]
    df <- as.integer(lavaan::fitmeasures(fit, "df"))
    expected <- ugamma_eigenvalues(
      lavaan::lavInspect(fit, "UGamma"), df
    )
    actual <- lavaan_lambdas(fit, df)$ug_biased
    expect_equal(actual, expected, tolerance = 1e-8, info = estimator)
    expect_true(
      categorical_valid_p(pvalues(fit, c("SB", "SS", "ALL", "PEBA4"))),
      info = estimator
    )
    info <- attr(pvalues(fit, "PEBA4"), "semtests")
    expect_equal(info$estimator_requested, estimator, info = estimator)
    expect_equal(info$spectrum_source, "lavaan UGamma", info = estimator)
  }
})

test_that("categorical SB and SS reproduce lavaan's public corrections", {
  variants <- list(
    DWLS = c("WLSM", "WLSMV"),
    ULS = c("ULSM", "ULSMV")
  )
  for (base in names(variants)) {
    fit_sb <- if (base == "ULS") {
      lavaan::cfa(
        categorical_test_model, categorical_fixture$data,
        ordered = categorical_fixture$ordered, estimator = "ULSM"
      )
    } else {
      categorical_base_fits[[variants[[base]][1L]]]
    }
    fit_ss <- categorical_base_fits[[variants[[base]][2L]]]
    expect_equal(
      as.numeric(pvalues(fit_sb, "SB")),
      lavaan::lavInspect(fit_sb, "test")$satorra.bentler$pvalue,
      tolerance = 1e-8,
      info = base
    )
    expect_equal(
      as.numeric(pvalues(fit_ss, "SS")),
      lavaan::lavInspect(fit_ss, "test")$scaled.shifted$pvalue,
      tolerance = 1e-8,
      info = base
    )
  }
})

test_that("categorical spectra cover groups, constraints, parameterizations, and mixed indicators", {
  cases <- list(
    multigroup = categorical_multigroup,
    constrained = categorical_pairs$default$m0,
    theta = categorical_pairs$theta$m1,
    mixed = categorical_pairs$mixed$m1,
    binary = {
      fixture <- categorical_test_data(binary = TRUE)
      lavaan::cfa(
        categorical_test_model, fixture$data,
        ordered = fixture$ordered
      )
    },
    pairwise = categorical_pairs$pairwise$m1
  )
  for (case in names(cases)) {
    fit <- cases[[case]]
    df <- as.integer(lavaan::fitmeasures(fit, "df"))
    expect_equal(
      lavaan_lambdas(fit, df)$ug_biased,
      ugamma_eigenvalues(lavaan::lavInspect(fit, "UGamma"), df),
      tolerance = 1e-8,
      info = case
    )
    expect_true(categorical_valid_p(pvalues(fit, c("SB", "SS", "PEBA4"))),
                info = case)
  }
})

test_that("full categorical WLS is refused (identity-spectrum no-op)", {
  fit <- fit_categorical_pair(estimator = "WLS")$m1
  df <- as.integer(lavaan::fitmeasures(fit, "df"))
  # The spectrum is still the identity -- that is exactly why it is refused ...
  expect_equal(lavaan_lambdas(fit, df)$ug_biased, rep(1, df), tolerance = 1e-8)
  # ... the public entry point declines rather than return the naive test.
  expect_error(pvalues(fit, c("STD", "SB", "SS", "ALL", "PEBA4")),
               "full weighted least squares")
})

test_that("nested categorical SB and SS reproduce lavaan Satorra-2000", {
  cases <- list(
    WLSMV = categorical_pairs$default,
    ULSMV = fit_categorical_pair("ULSMV"),
    pairwise = categorical_pairs$pairwise,
    mixed = categorical_pairs$mixed
  )
  for (case in names(cases)) {
    fits <- cases[[case]]
    scaled <- lavaan::lavTestLRT(
      fits$m1, fits$m0,
      method = "satorra.2000", A.method = "delta",
      scaled.shifted = FALSE
    )
    shifted <- lavaan::lavTestLRT(
      fits$m1, fits$m0,
      method = "satorra.2000", A.method = "delta",
      scaled.shifted = TRUE
    )
    actual <- pvalues_nested(
      fits$m0, fits$m1, tests = c("SB", "SS")
    )
    expect_equal(
      unname(actual[c("sb_ml", "ss_ml")]),
      c(
        unname(scaled[2L, "Pr(>Chisq)"]),
        unname(shifted[2L, "Pr(>Chisq)"])
      ),
      tolerance = 1e-8,
      info = case
    )
    expect_equal(attr(actual, "semtests")$A.method, "delta")
  }
})

test_that("nested multigroup categorical spectra reproduce lavaan traces", {
  h0 <- "
    visual  =~ x1 + c(a1, a2)*x2 + c(a1, a2)*x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  "
  m1 <- categorical_multigroup
  m0 <- lavaan::cfa(
    h0, categorical_fixture$data,
    ordered = categorical_fixture$ordered, group = "school"
  )
  actual <- pvalues_nested(m0, m1, tests = c("SB", "SS"))
  shifted <- lavaan::lavTestLRT(
    m1, m0, method = "satorra.2000",
    A.method = "delta", scaled.shifted = TRUE
  )
  # With WLSMV model tests, lavaan labels scaled.shifted = FALSE as a
  # mean-and-variance adjusted difference and changes its reference df. The
  # semTests SB uses the ordinary trace/df scaling. Fit the same DWLS models
  # through WLSM to obtain lavaan's public simple-scaled reference.
  m1_sb <- lavaan::cfa(
    categorical_test_model, categorical_fixture$data,
    ordered = categorical_fixture$ordered, group = "school", estimator = "WLSM"
  )
  m0_sb <- lavaan::cfa(
    h0, categorical_fixture$data,
    ordered = categorical_fixture$ordered, group = "school", estimator = "WLSM"
  )
  scaled <- lavaan::lavTestLRT(
    m1_sb, m0_sb, method = "satorra.2000",
    A.method = "delta", scaled.shifted = FALSE
  )
  actual_sb <- pvalues_nested(m0_sb, m1_sb, tests = "SB")
  expect_equal(as.numeric(actual_sb),
               unname(scaled[2L, "Pr(>Chisq)"]),
               tolerance = 1e-8)
  expect_equal(as.numeric(actual["ss_ml"]),
               unname(shifted[2L, "Pr(>Chisq)"]),
               tolerance = 1e-8)
})

test_that("nested categorical compatibility checks fail clearly", {
  fits <- categorical_pairs$default
  expect_error(pvalues_nested(fits$m0, fits$m1, method = "2001"),
               "method = \"2000\"")
  expect_error(pvalues_nested(fits$m0, fits$m1, A.method = "exact"),
               "A.method = \"delta\"")

  dwls <- fit_categorical_pair("DWLS")$m1
  expect_error(pvalues_nested(fits$m0, dwls), "same requested")

  theta <- categorical_pairs$theta$m1
  expect_error(pvalues_nested(fits$m0, theta), "same parameterization")

  different_data <- categorical_test_data()
  different_data$data$x1[1L] <- different_data$data$x1[2L]
  other <- lavaan::cfa(
    categorical_test_model, different_data$data,
    ordered = different_data$ordered
  )
  expect_error(pvalues_nested(fits$m0, other), "same raw data")

  continuous <- lavaan::cfa(
    categorical_test_model, lavaan::HolzingerSwineford1939,
    estimator = "MLM"
  )
  expect_error(pvalues_nested(fits$m0, continuous), "same data type")
})

test_that("categorical fit and pair compatibility guards cover the support surface", {
  fits <- categorical_pairs$default

  bad_estimator <- fits$m1
  bad_estimator@Options$estimator <- "ML"
  expect_error(check_supported(bad_estimator),
               "currently support lavaan's DWLS and ULS")

  bad_missing <- fits$m1
  bad_missing@Options$missing[1L] <- "ml"
  expect_error(check_supported(bad_missing),
               "listwise.*pairwise")

  pairwise <- categorical_pairs$pairwise$m1
  expect_error(check_categorical_nested_pair(fits$m0, pairwise),
               "same missing-data mode")

  mixed <- categorical_pairs$mixed$m1
  expect_error(check_categorical_nested_pair(fits$m0, mixed),
               "same observed and ordered variables")

  multigroup <- categorical_multigroup
  expect_error(check_categorical_nested_pair(fits$m0, multigroup),
               "same groups and group sample sizes")
})

test_that("categorical restriction rank failures are reported at the public entry point", {
  fits <- categorical_pairs$default
  testthat::local_mocked_bindings(
    get_a_matrix = function(...) matrix(0, 0L, 0L),
    .package = "semTests"
  )
  expect_error(
    pvalues_nested(fits$m0, fits$m1),
    "restriction rank.*does not match the df difference"
  )
})

test_that("UGamma eigenvalue validation distinguishes numerical instability", {
  expect_length(ugamma_eigenvalues(diag(2L), 0L), 0L)
  expect_error(
    ugamma_eigenvalues(matrix(letters[1:4], 2L), 2L),
    "square numeric"
  )
  expect_error(
    ugamma_eigenvalues(matrix(1:6, 2L), 2L),
    "square numeric"
  )
  expect_error(
    ugamma_eigenvalues(diag(2L), 3L),
    "dimension is smaller"
  )
  expect_error(
    ugamma_eigenvalues(diag(c(1, Inf)), 2L),
    "non-finite"
  )
  huge <- .Machine$double.xmax
  expect_error(
    ugamma_eigenvalues(matrix(c(huge, huge, huge, -huge), 2L), 2L),
    "non-finite eigenvalues"
  )
  expect_error(
    ugamma_eigenvalues(diag(c(2, -1e-10)), 2L),
    "rank"
  )
  expect_error(
    ugamma_eigenvalues(diag(c(2, -0.1)), 2L),
    "materially negative"
  )
  expect_error(
    ugamma_eigenvalues(matrix(c(0, -1, 1, 0), 2L), 2L),
    "complex"
  )
  expect_error(
    ugamma_eigenvalues(diag(c(1, 0)), 2L),
    "rank"
  )
  expect_error(
    lavaan_lambdas(NULL, 1L),
    "could not construct UGamma"
  )
})
