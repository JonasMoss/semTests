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

test_that("biased single-model spectra come directly from lavaan UGamma", {
  fixture <- categorical_test_data()
  for (estimator in c(
    "WLSMV", "DWLS", "WLSM", "WLSMVS", "ULSMV", "ULS", "WLS"
  )) {
    fit <- lavaan::cfa(
      categorical_test_model, fixture$data,
      ordered = fixture$ordered, estimator = estimator
    )
    df <- as.integer(lavaan::fitmeasures(fit, "df"))
    expected <- ugamma_eigenvalues(
      lavaan::lavInspect(fit, "UGamma"), df
    )
    actual <- lavaan_lambdas(fit, df)$ug_biased
    expect_equal(actual, expected, tolerance = 1e-12, info = estimator)
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
  fixture <- categorical_test_data()
  variants <- list(
    DWLS = c("WLSM", "WLSMV"),
    ULS = c("ULSM", "ULSMV")
  )
  for (base in names(variants)) {
    fit_sb <- lavaan::cfa(
      categorical_test_model, fixture$data,
      ordered = fixture$ordered, estimator = variants[[base]][1L]
    )
    fit_ss <- lavaan::cfa(
      categorical_test_model, fixture$data,
      ordered = fixture$ordered, estimator = variants[[base]][2L]
    )
    expect_equal(
      as.numeric(pvalues(fit_sb, "SB")),
      lavaan::lavInspect(fit_sb, "test")$satorra.bentler$pvalue,
      tolerance = 1e-12,
      info = base
    )
    expect_equal(
      as.numeric(pvalues(fit_ss, "SS")),
      lavaan::lavInspect(fit_ss, "test")$scaled.shifted$pvalue,
      tolerance = 1e-12,
      info = base
    )
  }
})

test_that("categorical spectra cover groups, constraints, parameterizations, and mixed indicators", {
  multigroup_fixture <- categorical_test_data()
  cases <- list(
    multigroup = lavaan::cfa(
      categorical_test_model, multigroup_fixture$data,
      ordered = multigroup_fixture$ordered, group = "school"
    ),
    constrained = fit_categorical_pair()$m0,
    theta = fit_categorical_pair(parameterization = "theta")$m1,
    mixed = fit_categorical_pair(mixed = TRUE)$m1,
    binary = {
      fixture <- categorical_test_data(binary = TRUE)
      lavaan::cfa(
        categorical_test_model, fixture$data,
        ordered = fixture$ordered
      )
    },
    pairwise = fit_categorical_pair(pairwise = TRUE)$m1
  )
  for (case in names(cases)) {
    fit <- cases[[case]]
    df <- as.integer(lavaan::fitmeasures(fit, "df"))
    expect_equal(
      lavaan_lambdas(fit, df)$ug_biased,
      ugamma_eigenvalues(lavaan::lavInspect(fit, "UGamma"), df),
      tolerance = 1e-12,
      info = case
    )
    expect_true(categorical_valid_p(pvalues(fit, c("SB", "SS", "PEBA4"))),
                info = case)
  }
})

test_that("full categorical WLS is the identity-spectrum no-op", {
  fit <- fit_categorical_pair(estimator = "WLS")$m1
  df <- as.integer(lavaan::fitmeasures(fit, "df"))
  spectrum <- lavaan_lambdas(fit, df)$ug_biased
  reference <- as.numeric(
    1 - stats::pchisq(lavaan::fitmeasures(fit, "chisq"), df)
  )
  expect_equal(spectrum, rep(1, df), tolerance = 1e-8)
  expect_equal(
    as.numeric(pvalues(fit, c("STD", "SB", "SS", "ALL", "PEBA4"))),
    rep(reference, 5L),
    tolerance = 1e-8
  )
})

test_that("nested categorical SB and SS reproduce lavaan Satorra-2000", {
  cases <- list(
    WLSMV = fit_categorical_pair("WLSMV"),
    ULSMV = fit_categorical_pair("ULSMV"),
    pairwise = fit_categorical_pair("WLSMV", pairwise = TRUE),
    mixed = fit_categorical_pair("WLSMV", mixed = TRUE)
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
      tolerance = 1e-10,
      info = case
    )
    expect_equal(attr(actual, "semtests")$A.method, "delta")
  }
})

test_that("nested multigroup categorical spectra reproduce lavaan traces", {
  fixture <- categorical_test_data()
  h0 <- "
    visual  =~ x1 + c(a1, a2)*x2 + c(a1, a2)*x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  "
  m1 <- lavaan::cfa(
    categorical_test_model, fixture$data,
    ordered = fixture$ordered, group = "school"
  )
  m0 <- lavaan::cfa(
    h0, fixture$data,
    ordered = fixture$ordered, group = "school"
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
    categorical_test_model, fixture$data,
    ordered = fixture$ordered, group = "school", estimator = "WLSM"
  )
  m0_sb <- lavaan::cfa(
    h0, fixture$data,
    ordered = fixture$ordered, group = "school", estimator = "WLSM"
  )
  scaled <- lavaan::lavTestLRT(
    m1_sb, m0_sb, method = "satorra.2000",
    A.method = "delta", scaled.shifted = FALSE
  )
  actual_sb <- pvalues_nested(m0_sb, m1_sb, tests = "SB")
  expect_equal(as.numeric(actual_sb),
               unname(scaled[2L, "Pr(>Chisq)"]),
               tolerance = 1e-10)
  expect_equal(as.numeric(actual["ss_ml"]),
               unname(shifted[2L, "Pr(>Chisq)"]),
               tolerance = 1e-10)
})

test_that("nested categorical compatibility checks fail clearly", {
  fits <- fit_categorical_pair()
  expect_error(pvalues_nested(fits$m0, fits$m1, method = "2001"),
               "method = \"2000\"")
  expect_error(pvalues_nested(fits$m0, fits$m1, A.method = "exact"),
               "A.method = \"delta\"")

  dwls <- fit_categorical_pair("DWLS")$m1
  expect_error(pvalues_nested(fits$m0, dwls), "same requested")

  theta <- fit_categorical_pair(parameterization = "theta")$m1
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
  fits <- fit_categorical_pair()

  bad_estimator <- fits$m1
  bad_estimator@Options$estimator <- "ML"
  expect_error(check_supported(bad_estimator),
               "currently support lavaan's DWLS, ULS, and WLS")

  bad_missing <- fits$m1
  bad_missing@Options$missing[1L] <- "ml"
  expect_error(check_supported(bad_missing),
               "listwise.*pairwise")

  pairwise <- fit_categorical_pair(pairwise = TRUE)$m1
  expect_error(check_categorical_nested_pair(fits$m0, pairwise),
               "same missing-data mode")

  fixture <- categorical_test_data(mixed = TRUE)
  mixed <- lavaan::cfa(
    categorical_test_model, fixture$data,
    ordered = fixture$ordered
  )
  expect_error(check_categorical_nested_pair(fits$m0, mixed),
               "same observed and ordered variables")

  fixture <- categorical_test_data()
  multigroup <- lavaan::cfa(
    categorical_test_model, fixture$data,
    ordered = fixture$ordered, group = "school"
  )
  expect_error(check_categorical_nested_pair(fits$m0, multigroup),
               "same groups and group sample sizes")
})

test_that("categorical restriction rank failures are reported at the public entry point", {
  fits <- fit_categorical_pair()
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
