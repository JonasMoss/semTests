# The eigenvalue-based p-values are valid for any minimum-discrepancy estimator
# (ADF excepted, where they reduce to the exact chi-square). These tests
# exercise the broadened estimator / data-type support: single-model for all
# estimators, nested for continuous estimators, and the guards that keep the
# normal-theory-only pieces (RLS statistic, Du-Bentler unbiased gamma) and
# unsupported categorical nested variants from being used.

hs  <- " visual =~ x1 + x2 + x3
         textual =~ x4 + x5 + x6
         speed =~ x7 + x8 + x9 "
hs0 <- " visual =~ x1 + a*x2 + x3
         textual =~ x4 + a*x5 + a*x6
         speed =~ x7 + a*x8 + x9 "

ordinalize <- function(data, vars = paste0("x", 1:9)) {
  for (v in vars) data[[v]] <- ordered(cut(data[[v]], 3))
  data
}

valid_p <- function(p) all(is.finite(p) & p >= 0 & p <= 1)

test_that("single-model p-values work for ULS / GLS / MLR / FIML", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS
  set.seed(1)
  HSm$x1[sample(nrow(HS), 60)] <- NA
  fits <- list(
    ULS  = lavaan::cfa(hs, HS,  estimator = "ULS", test = "satorra.bentler"),
    GLS  = lavaan::cfa(hs, HS,  estimator = "GLS"),
    MLR  = lavaan::cfa(hs, HS,  estimator = "MLR"),
    FIML = lavaan::cfa(hs, HSm, missing = "fiml", estimator = "MLR")
  )
  for (nm in names(fits)) {
    p <- pvalues(fits[[nm]], "PEBA4")
    expect_true(valid_p(p), info = nm)
    expect_s3_class(p, "semTests_pvalues")
    expect_equal(attr(p, "semtests")$estimator, fits[[nm]]@Options$estimator, info = nm)
  }
})

test_that("single-model p-values work for categorical WLSMV", {
  fit <- lavaan::cfa(hs, ordinalize(lavaan::HolzingerSwineford1939),
                     ordered = paste0("x", 1:9))
  p <- pvalues(fit, "PEBA4")
  expect_true(valid_p(p))
  expect_equal(attr(p, "semtests")$data_type, "categorical")
})

test_that("nested p-values work for continuous non-ML estimators", {
  HS <- lavaan::HolzingerSwineford1939
  # ULS needs a robust test for the NACOV (gamma); GLS supplies it natively and
  # rejects test = "satorra.bentler".
  uls1 <- lavaan::cfa(hs,  HS, estimator = "ULS", test = "satorra.bentler")
  uls0 <- lavaan::cfa(hs0, HS, estimator = "ULS", test = "satorra.bentler")
  expect_true(valid_p(pvalues_nested(uls0, uls1)), info = "ULS")
  gls1 <- lavaan::cfa(hs,  HS, estimator = "GLS")
  gls0 <- lavaan::cfa(hs0, HS, estimator = "GLS")
  expect_true(valid_p(pvalues_nested(gls0, gls1)), info = "GLS")
})

test_that("normal-theory-only options are refused off the normal-theory ML case", {
  # is_classic_nt() keys on the *estimator* (ML), so MLM/MLR on complete data
  # still admit RLS/UG; the non-classic fits are GLS, ULS, categorical, and FIML.
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(8); HSm$x1[sample(nrow(HS), 40)] <- NA
  non_classic <- list(
    ULS         = lavaan::cfa(hs, HS, estimator = "ULS", test = "satorra.bentler"),
    GLS         = lavaan::cfa(hs, HS, estimator = "GLS"),
    categorical = lavaan::cfa(hs, ordinalize(HS), ordered = paste0("x", 1:9)),
    FIML        = lavaan::cfa(hs, HSm, missing = "fiml", estimator = "MLR")
  )
  for (nm in names(non_classic)) {
    expect_error(pvalues(non_classic[[nm]], "SB_RLS"), "RLS",      info = nm)
    expect_error(pvalues(non_classic[[nm]], "SB_UG"),  "unbiased", info = nm)
  }
})

test_that("FIML constraints (single-group, no fixed exogenous covariates) are enforced", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(9); HSm$x1[sample(nrow(HS), 40)] <- NA
  mg_fiml <- lavaan::cfa(hs, HSm, missing = "fiml", estimator = "MLR",
                         group = "school")
  expect_error(pvalues(mg_fiml, "PEBA4"), "single-group")

  fixed_x <- lavaan::cfa(
    "visual =~ x1 + x2 + x3
     visual ~ ageyr",
    HSm, missing = "fiml", estimator = "MLR", fixed.x = TRUE
  )
  expect_error(pvalues(fixed_x, "PEBA4"), "fixed exogenous")
})

test_that("nested categorical works through the 2000/delta path", {
  HSc <- ordinalize(lavaan::HolzingerSwineford1939)
  fc1 <- lavaan::cfa(hs,  HSc, ordered = paste0("x", 1:9))
  fc0 <- lavaan::cfa(hs0, HSc, ordered = paste0("x", 1:9))
  expect_true(valid_p(pvalues_nested(fc0, fc1)))
  expect_error(pvalues_nested(fc0, fc1, method = "2001"),
               "method = \"2000\"")
  expect_error(pvalues_nested(fc0, fc1, A.method = "exact"),
               "A.method = \"delta\"")
})

test_that("nested missing-data (FIML) works for continuous fits", {
  HS <- lavaan::HolzingerSwineford1939
  set.seed(2); HS$x1[sample(nrow(HS), 40)] <- NA
  f1 <- lavaan::cfa(hs,  HS, missing = "fiml", estimator = "MLR")
  f0 <- lavaan::cfa(hs0, HS, missing = "fiml", estimator = "MLR")
  p_delta <- pvalues_nested(f0, f1)
  p_exact <- pvalues_nested(f0, f1, A.method = "exact")
  expect_true(valid_p(p_exact))
  expect_true(valid_p(p_delta))
  expect_equal(attr(p_delta, "semtests")$A.method, "delta")
  expect_equal(attr(p_exact, "semtests")$A.method, "exact")
  expect_error(pvalues_nested(f0, f1, method = "2001"), "method = \"2000\"")
  expect_error(pvalues_nested(f0, f1, tests = "PALL_UG"), "unbiased")
})

test_that("the output object records the options it used", {
  fit  <- lavaan::cfa(hs, lavaan::HolzingerSwineford1939, estimator = "MLM")
  info <- attr(pvalues(fit, c("PEBA4_RLS", "SB_UG_ML")), "semtests")
  expect_equal(info$estimator, "ML")
  expect_equal(info$data_type, "continuous")
  expect_false(info$nested)
  expect_true(is.numeric(info$df))
  expect_equal(info$requested_tests, c("PEBA4_RLS", "SB_UG_ML"))
  expect_equal(unname(info$statistic), c("rls", "ml"))
  expect_equal(unname(info$gamma), c("biased", "unbiased"))
})

test_that("non-lavaan objects are rejected with a clear message", {
  fit <- lavaan::cfa(hs, lavaan::HolzingerSwineford1939, estimator = "MLM")
  # Single-model: NULL, a list, and a data.frame must all be refused up front.
  expect_error(pvalues(NULL), "lavaan")
  expect_error(pvalues(list(a = 1)), "lavaan")
  expect_error(pvalues(lavaan::HolzingerSwineford1939), "lavaan")
  # Nested: a non-lavaan in either slot must fail by name *before* the df
  # computation reads m0@test / m1@test (the ordering regression).
  expect_error(pvalues_nested(NULL, fit), "m0")
  expect_error(pvalues_nested(fit, NULL), "m1")
  expect_error(pvalues_nested(fit, "not a fit"), "m1")
})
