fiml_test_data <- function() {
  data <- lavaan::HolzingerSwineford1939
  for (j in seq_len(9L)) {
    data[
      seq(3L + j, nrow(data), by = 17L + j),
      paste0("x", j)
    ] <- NA_real_
  }
  data
}

fiml_h1_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

fiml_h0_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + a*x5 + a*x6
  speed   =~ x7 + x8 + x9
"

fit_fiml_pair <- function() {
  data <- fiml_test_data()
  list(
    m0 = lavaan::cfa(
      fiml_h0_model, data, missing = "fiml", estimator = "MLR",
      meanstructure = TRUE
    ),
    m1 = lavaan::cfa(
      fiml_h1_model, data, missing = "fiml", estimator = "MLR",
      meanstructure = TRUE
    )
  )
}

test_that("lavaan FIML single-model convention is its inspected UGamma spectrum", {
  fit <- fit_fiml_pair()$m1
  df <- as.integer(lavaan::fitmeasures(fit, "df"))
  actual <- fiml_lambdas(fit, df, fiml.convention = "lavaan")$ug_biased
  expected <- sort(
    Re(eigen(lavaan::lavInspect(fit, "UGamma"),
             only.values = TRUE)$values),
    decreasing = TRUE
  )[seq_len(df)]
  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("lavaan nested convention reproduces public Satorra-2000 tests", {
  fits <- fit_fiml_pair()
  for (A.method in c("delta", "exact")) {
    scaled <- lavaan::lavTestLRT(
      fits$m1, fits$m0,
      method = "satorra.2000",
      A.method = A.method,
      scaled.shifted = FALSE
    )
    shifted <- lavaan::lavTestLRT(
      fits$m1, fits$m0,
      method = "satorra.2000",
      A.method = A.method,
      scaled.shifted = TRUE
    )
    actual <- pvalues_nested(
      fits$m0, fits$m1,
      tests = c("SB", "SS"),
      A.method = A.method,
      fiml.convention = "lavaan"
    )
    expect_equal(
      unname(actual[c("sb_ml", "ss_ml")]),
      c(
        unname(scaled[2L, "Pr(>Chisq)"]),
        unname(shifted[2L, "Pr(>Chisq)"])
      ),
      tolerance = 1e-8,
      info = A.method
    )
  }
})

test_that("observed FIML spectra match pinned magmaan validation goldens", {
  # Generated with magmaan commit
  # ceacbf5b7ab35353d1ade0cf956ff7d0dff256d5. magmaan is deliberately not a
  # runtime/test dependency: these fixed values make the independent
  # implementation an oracle without making CRAN install another package.
  expected_h1 <- c(
    0.516997470209856, 0.533197489524983, 0.612566133388422,
    0.649910769397448, 0.658426886876557, 0.686454389363300,
    0.727518146307215, 0.742999997015816, 0.805203906207252,
    0.858190223354869, 0.926566717574974, 0.950822887211618,
    1.028363016942060, 1.098471402533230, 1.113892184586230,
    1.180704681887030, 1.240329629090810, 1.277264374290930,
    1.356741223677690, 1.411144905844780, 1.440443613049470,
    1.526601595880390, 1.673403765187960, 1.788877181382700
  )
  expected_h0 <- c(
    0.493699915318542, 0.519247884618317, 0.589996912538212,
    0.625386027537290, 0.658802419650294, 0.683286346347235,
    0.727447426051526, 0.756567246294660, 0.783765077875271,
    0.857446878975972, 0.878857510678687, 0.950873268356167,
    1.032127336741360, 1.084902779011540, 1.106458183138980,
    1.154264988900710, 1.193672538436130, 1.262692787447020,
    1.349871075104820, 1.396636138702200, 1.427960209477430,
    1.479304541220980, 1.539643425180680, 1.658317648201370,
    1.813988710938970
  )
  fits <- fit_fiml_pair()
  h1 <- sort(fiml_lambdas(
    fits$m1, 24L, fiml.convention = "observed"
  )$ug_biased)
  h0 <- sort(fiml_lambdas(
    fits$m0, 25L, fiml.convention = "observed"
  )$ug_biased)
  delta <- fiml_lambdas_nested(
    fits$m0, fits$m1, 1L,
    A.method = "delta", fiml.convention = "observed"
  )$ug_biased
  exact <- fiml_lambdas_nested(
    fits$m0, fits$m1, 1L,
    A.method = "exact", fiml.convention = "observed"
  )$ug_biased

  expect_equal(h1, expected_h1, tolerance = 1e-4)
  expect_equal(h0, expected_h0, tolerance = 1e-4)
  expect_equal(delta, 1.36384308952366, tolerance = 1e-4)
  expect_equal(exact, 1.44378023818416, tolerance = 1e-4)
})

test_that("FIML equality bases cover labels, general constraints, and mixtures", {
  data <- fiml_test_data()
  shared <- "
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  "
  models <- list(
    simple = paste("visual =~ x1 + a*x2 + a*x3", shared),
    general = paste(
      "visual =~ x1 + a*x2 + b*x3", shared, "a == b",
      sep = "\n"
    ),
    mixed = paste(
      "visual =~ x1 + a*x2 + b*x3", shared, "a == b", "a > .1",
      sep = "\n"
    )
  )
  spectra <- lapply(models, function(model) {
    fit <- lavaan::cfa(model, data, missing = "fiml", estimator = "MLR")
    K <- fiml_K_matrix(fit)
    expect_equal(nrow(K), fit@Model@nx.free)
    expect_equal(ncol(K), fit@Model@nx.free - 1L)
    expect_equal(crossprod(K), diag(ncol(K)), tolerance = 1e-8)
    fiml_lambdas(
      fit, as.integer(lavaan::fitmeasures(fit, "df")), "observed"
    )$ug_biased
  })
  expect_equal(spectra$simple, spectra$general, tolerance = 1e-8)
  expect_equal(spectra$simple, spectra$mixed, tolerance = 3e-5)
})

test_that("nested FIML supports an already-constrained H1", {
  data <- fiml_test_data()
  h1 <- "
    visual  =~ x1 + a*x2 + b*x3
    textual =~ x4 + c*x5 + d*x6
    speed   =~ x7 + x8 + x9
    a == b
  "
  h0 <- paste(h1, "c == d", sep = "\n")
  m1 <- lavaan::cfa(h1, data, missing = "fiml", estimator = "MLR")
  m0 <- lavaan::cfa(h0, data, missing = "fiml", estimator = "MLR")
  expect_equal(ncol(fiml_K_matrix(m1)), m1@Model@nx.free - 1L)
  for (A.method in c("delta", "exact")) {
    spectrum <- fiml_lambdas_nested(
      m0, m1, 1L, A.method = A.method,
      fiml.convention = "observed"
    )$ug_biased
    expect_length(spectrum, 1L)
    expect_true(is.finite(spectrum) && spectrum > 0)
  }
})

test_that("nested FIML restrictions are invariant to parameter-table ordering", {
  data <- fiml_test_data()
  h1 <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
    visual ~~ textual
    visual ~~ speed
    textual ~~ speed
  "
  h0_aligned <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + a*x5 + a*x6
    speed   =~ x7 + x8 + x9
    visual ~~ textual
    visual ~~ speed
    textual ~~ speed
  "
  h0_reordered <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + a*x5 + a*x6
    speed   =~ x7 + x8 + x9
    textual ~~ speed
    visual ~~ speed
    visual ~~ textual
  "
  m1 <- lavaan::cfa(h1, data, missing = "fiml", estimator = "MLR")
  m0a <- lavaan::cfa(h0_aligned, data, missing = "fiml", estimator = "MLR")
  m0b <- lavaan::cfa(h0_reordered, data, missing = "fiml", estimator = "MLR")
  expect_false(identical(fiml_parameter_keys(m0a), fiml_parameter_keys(m0b)))
  for (A.method in c("delta", "exact")) {
    a <- fiml_lambdas_nested(m0a, m1, 1L, A.method, "observed")$ug_biased
    b <- fiml_lambdas_nested(m0b, m1, 1L, A.method, "observed")$ug_biased
    expect_equal(a, b, tolerance = 1e-8, info = A.method)
  }
})

test_that("FIML delta supports parameter-set changes that exact rejects clearly", {
  data <- fiml_test_data()
  h1 <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  "
  h0 <- paste(h1, "visual ~~ 0*textual", sep = "\n")
  m1 <- lavaan::cfa(h1, data, missing = "fiml", estimator = "MLR")
  m0 <- lavaan::cfa(h0, data, missing = "fiml", estimator = "MLR")

  delta <- fiml_lambdas_nested(
    m0, m1, 1L,
    A.method = "delta", fiml.convention = "observed"
  )$ug_biased
  expect_length(delta, 1L)
  expect_true(is.finite(delta) && delta > 0)
  expect_error(
    fiml_lambdas_nested(
      m0, m1, 1L,
      A.method = "exact", fiml.convention = "observed"
    ),
    "same underlying full free-parameter set"
  )
})
