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
      fiml_h0_model, data,
      missing = "fiml", estimator = "MLR",
      meanstructure = TRUE
    ),
    m1 = lavaan::cfa(
      fiml_h1_model, data,
      missing = "fiml", estimator = "MLR",
      meanstructure = TRUE
    )
  )
}

fit_fiml_multigroup_pair <- function() {
  data <- fiml_test_data()
  list(
    m0 = lavaan::cfa(
      fiml_h1_model, data,
      group = "school", group.equal = "loadings",
      missing = "fiml", estimator = "MLR", meanstructure = TRUE
    ),
    m1 = lavaan::cfa(
      fiml_h1_model, data,
      group = "school",
      missing = "fiml", estimator = "MLR", meanstructure = TRUE
    )
  )
}

fit_fiml_random_x_pair <- function() {
  data <- fiml_test_data()
  data$ageyr[seq(7L, nrow(data), by = 23L)] <- NA_real_
  h1 <- "
    visual =~ x1 + x2 + x3
    visual ~ ageyr
  "
  h0 <- "
    visual =~ x1 + a*x2 + a*x3
    visual ~ ageyr
  "
  list(
    m0 = lavaan::cfa(
      h0, data,
      missing = "fiml", estimator = "MLR",
      fixed.x = FALSE, conditional.x = FALSE, meanstructure = TRUE
    ),
    m1 = lavaan::cfa(
      h1, data,
      missing = "fiml", estimator = "MLR",
      fixed.x = FALSE, conditional.x = FALSE, meanstructure = TRUE
    )
  )
}

fiml_pair <- fit_fiml_pair()
fiml_multigroup_pair <- fit_fiml_multigroup_pair()
fiml_random_x_pair <- fit_fiml_random_x_pair()

test_that("lavaan FIML single-model convention is its inspected UGamma spectrum", {
  fit <- fiml_pair$m1
  df <- as.integer(lavaan::fitmeasures(fit, "df"))
  actual <- fiml_lambdas(fit, df, fiml.convention = "lavaan")$ug_biased
  expected <- sort(
    Re(eigen(lavaan::lavInspect(fit, "UGamma"),
      only.values = TRUE
    )$values),
    decreasing = TRUE
  )[seq_len(df)]
  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("lavaan nested convention reproduces public Satorra-2000 tests", {
  fits <- fiml_pair
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

test_that("multigroup lavaan convention reproduces inspected and public tests", {
  fits <- fiml_multigroup_pair
  df <- as.integer(lavaan::fitmeasures(fits$m1, "df"))
  actual <- fiml_lambdas(
    fits$m1, df,
    fiml.convention = "lavaan"
  )$ug_biased
  expected <- sort(
    Re(eigen(lavaan::lavInspect(fits$m1, "UGamma"),
      only.values = TRUE
    )$values),
    decreasing = TRUE
  )[seq_len(df)]
  expect_equal(actual, expected, tolerance = 1e-8)

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
    semtests <- pvalues_nested(
      fits$m0, fits$m1,
      tests = c("SB", "SS"),
      A.method = A.method,
      fiml.convention = "lavaan"
    )
    expect_equal(
      unname(semtests[c("sb_ml", "ss_ml")]),
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
  fits <- fiml_pair
  h1 <- sort(fiml_lambdas(
    fits$m1, 24L,
    fiml.convention = "observed"
  )$ug_biased)
  h0 <- sort(fiml_lambdas(
    fits$m0, 25L,
    fiml.convention = "observed"
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

test_that("multigroup observed FIML matches pinned magmaan goldens", {
  # Generated with magmaan commit
  # a48a177cd50466c79e5a214beea5c0297386881d. The two school groups have
  # unequal sizes, and the nested pair imposes cross-group loading equalities.
  expected_single <- c(
    2.348991346666434, 2.144946798341949, 2.125122259230327,
    1.896750302736278, 1.895093422434539, 1.732796134720916,
    1.667194253006357, 1.639963324912957, 1.611119755448654,
    1.445125781254547, 1.399273927084085, 1.345691892141516,
    1.309931846918966, 1.279272550851228, 1.138342023694073,
    1.119628477703533, 1.068729115127381, 1.041070100684627,
    1.026516221652800, 1.023382139605715, 0.976765876188130,
    0.957381912612311, 0.926709286231006, 0.906882388301446,
    0.896751425741604, 0.840998953825102, 0.826759425794660,
    0.795725471575969, 0.783635287333617, 0.752746278156110,
    0.708242746447536, 0.686945213162457, 0.668667830229329,
    0.595607184795556, 0.585803522770502, 0.565619965235117,
    0.560577245332923, 0.533518442137667, 0.523044659059982,
    0.509485764364486, 0.479503542489589, 0.470272174367425,
    0.432302072557310, 0.415950165250539, 0.382940197188285,
    0.380028745045868, 0.340733089387116, 0.301773187780947
  )
  expected_delta <- c(
    1.869804285582957, 1.702368648168416, 1.640880437868678,
    1.207303228193902, 1.010229424023672, 0.887385866841645
  )
  expected_exact <- c(
    1.983097677719228, 1.767310942211078, 1.712195417641462,
    1.149782256544814, 0.998486667107969, 0.876532808180277
  )

  fits <- fiml_multigroup_pair
  df <- as.integer(lavaan::fitmeasures(fits$m1, "df"))
  single <- fiml_lambdas(
    fits$m1, df,
    fiml.convention = "observed"
  )$ug_biased
  nested_df <- as.integer(
    lavaan::fitmeasures(fits$m0, "df") -
      lavaan::fitmeasures(fits$m1, "df")
  )
  delta <- fiml_lambdas_nested(
    fits$m0, fits$m1, nested_df,
    A.method = "delta", fiml.convention = "observed"
  )$ug_biased
  exact <- fiml_lambdas_nested(
    fits$m0, fits$m1, nested_df,
    A.method = "exact", fiml.convention = "observed"
  )$ug_biased

  expect_equal(single, expected_single, tolerance = 1e-4)
  expect_equal(delta, expected_delta, tolerance = 3e-4)
  expect_equal(exact, expected_exact, tolerance = 3e-4)
  expect_true(all(is.finite(pvalues(fits$m1, c("SB", "SS", "PEBA4")))))
  expect_true(all(is.finite(pvalues_nested(
    fits$m0, fits$m1,
    tests = c("SB", "SS", "PEBA4")
  ))))
})

test_that("random-x FIML matches pinned magmaan goldens", {
  # Generated with magmaan commit
  # a48a177cd50466c79e5a214beea5c0297386881d. The observed predictor is
  # modeled jointly by setting fixed.x = FALSE in both implementations.
  expected_single <- c(
    1.143360126227509, 1.015000221043247
  )
  expected_delta <- 0.906390628017868
  expected_exact <- 1.060479849465168

  fits <- fiml_random_x_pair
  df <- as.integer(lavaan::fitmeasures(fits$m1, "df"))
  nested_df <- as.integer(
    lavaan::fitmeasures(fits$m0, "df") -
      lavaan::fitmeasures(fits$m1, "df")
  )
  single <- fiml_lambdas(
    fits$m1, df,
    fiml.convention = "observed"
  )$ug_biased
  delta <- fiml_lambdas_nested(
    fits$m0, fits$m1, nested_df,
    A.method = "delta", fiml.convention = "observed"
  )$ug_biased
  exact <- fiml_lambdas_nested(
    fits$m0, fits$m1, nested_df,
    A.method = "exact", fiml.convention = "observed"
  )$ug_biased

  expect_equal(single, expected_single, tolerance = 1e-5)
  expect_equal(delta, expected_delta, tolerance = 1e-5)
  expect_equal(exact, expected_exact, tolerance = 1e-5)
  expect_true(all(is.finite(pvalues(fits$m1, c("SB", "SS", "PEBA2")))))
})

test_that("random-x FIML lavaan convention reproduces public nested tests", {
  fits <- fiml_random_x_pair
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
    K <- model_parameter_basis(fit)
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
  expect_equal(ncol(model_parameter_basis(m1)), m1@Model@nx.free - 1L)
  for (A.method in c("delta", "exact")) {
    spectrum <- fiml_lambdas_nested(
      m0, m1, 1L,
      A.method = A.method,
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
