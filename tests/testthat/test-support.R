# Conformance tests for the documented support matrix (?semTests-support) and
# the entry-level gate (check_supported / check_supported_nested). The matrix is
# the source of truth: every available row should compute, every rejected row
# should stop, and each FIML information convention should be explicit through
# the public entry point. setup.R supplies the multi-group fixtures (m0/m1 MLM,
# m0_/m1_ GLS, object, m1_no_groups/m0_no_groups).

hs <- " visual =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6
        speed =~ x7 + x8 + x9 "

ordinalize <- function(data, vars = paste0("x", 1:9)) {
  for (v in vars) data[[v]] <- ordered(cut(data[[v]], 3))
  data
}
valid_p <- function(p) all(is.finite(p) & p >= 0 & p <= 1)

test_that("multi-group continuous single-model p-values are supported", {
  p_mlm <- pvalues(m1, "PEBA4")       # MLM, group = "school"
  p_gls <- pvalues(m1_, "PEBA4")      # GLS, group = "school"
  expect_true(valid_p(p_mlm))
  expect_true(valid_p(p_gls))
  expect_equal(attr(p_mlm, "semtests")$estimator, "ML")
  expect_equal(attr(p_gls, "semtests")$estimator, "GLS")
  expect_equal(attr(p_mlm, "semtests")$data_type, "continuous")
})

test_that("multi-group continuous nested p-values are supported (both methods)", {
  expect_true(valid_p(pvalues_nested(m0,  m1)))                    # MLM, 2000
  expect_true(valid_p(pvalues_nested(m0,  m1, method = "2001")))   # MLM, 2001
  expect_true(valid_p(pvalues_nested(m0_, m1_)))                   # GLS, 2000
})

test_that("multi-group categorical single-model is supported", {
  mgc <- lavaan::cfa(hs, ordinalize(lavaan::HolzingerSwineford1939),
                     ordered = paste0("x", 1:9), group = "school")
  p <- pvalues(mgc, "PEBA4")
  expect_true(valid_p(p))
  expect_equal(attr(p, "semtests")$data_type, "categorical")
})

test_that("ADF/WLS is degenerate: every p-value equals the ordinary chi-square", {
  skip_on_cran()
  set.seed(2024)
  dat <- lavaan::simulateData(hs, sample.nobs = 600)
  wls <- lavaan::cfa(hs, dat, estimator = "WLS")
  ref <- as.numeric(1 - stats::pchisq(lavaan::fitmeasures(wls, "chisq"),
                                       lavaan::fitmeasures(wls, "df")))
  p <- pvalues(wls, c("STD", "PEBA4", "SB", "POLS"))
  expect_true(valid_p(p))
  expect_equal(unname(as.numeric(p)), rep(ref, length(p)), tolerance = 1e-8)
})

test_that("check_supported rejects estimators outside the supported set", {
  bad <- object
  bad@Options$estimator <- "PML"
  expect_error(check_supported(bad), "not supported")
})

test_that("check_supported rejects non-FIML missing data", {
  bad <- object
  bad@Options$missing <- "pairwise"
  expect_error(check_supported(bad), "FIML")
})

test_that("nested missing-data requires FIML for both fits", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(11); HSm$x1[sample(nrow(HS), 40)] <- NA
  fiml     <- lavaan::cfa(hs, HSm, missing = "fiml", estimator = "MLR")
  complete <- lavaan::cfa(hs, HS,  estimator = "MLM")
  expect_error(check_supported_nested(fiml, complete, "2000", "exact"), "FIML")
})

test_that("the FIML convention is explicit and recorded", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(12); HSm$x1[sample(nrow(HS), 60)] <- NA
  fiml_expected <- lavaan::cfa(hs, HSm, missing = "fiml", estimator = "ML",
                               information = "expected")

  observed <- pvalues(fiml_expected, "PEBA4")
  lavaan <- pvalues(fiml_expected, "PEBA4", fiml.convention = "lavaan")
  expect_equal(attr(observed, "semtests")$fiml.convention, "observed")
  expect_equal(attr(lavaan, "semtests")$fiml.convention, "lavaan")
  expect_null(attr(pvalues(object, "PEBA4_RLS"), "semtests")$fiml.convention)
})
