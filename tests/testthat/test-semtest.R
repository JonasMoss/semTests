# semtest() records the settings that produced each result, under the names a
# user passes, so a saved result says how it was computed.

semtest_model <- "f =~ x1 + a*x2 + b*x3 + c*x4"
semtest_fit <- lavaan::cfa(semtest_model, lavaan::HolzingerSwineford1939,
                           meanstructure = TRUE)
semtest_restricted <- lavaan::cfa(
  paste(semtest_model, "a == b", "b == c", sep = "\n"),
  lavaan::HolzingerSwineford1939, meanstructure = TRUE
)

test_that("nested results record the resolved settings as labels", {
  out <- semtest(semtest_restricted, semtest_fit)
  expect_identical(out$provenance$curvature, "observed")
  expect_identical(out$provenance$score_covariance, "centered")
  expect_identical(out$provenance$weight, "curvature")

  out <- semtest(semtest_restricted, semtest_fit, weight = "expected_curvature",
                 curvature = "expected", score_covariance = "uncentered")
  expect_identical(out$provenance$curvature, "expected")
  expect_identical(out$provenance$score_covariance, "uncentered")
  expect_identical(out$provenance$weight, "expected_curvature")

  out <- semtest(semtest_restricted, semtest_fit, statistic = "criterion",
                 weight = "score_covariance")
  expect_identical(out$type, "criterion difference")
  expect_identical(out$provenance$curvature, "observed")
  expect_identical(out$provenance$weight, "curvature")
})

test_that("overall-fit results record the resolved settings as labels", {
  out <- semtest(semtest_fit)
  expect_identical(out$provenance$curvature, "expected")
  expect_identical(out$provenance$score_covariance, "uncentered")
  expect_identical(out$provenance$weight, "curvature")
  out <- semtest(semtest_fit, weight = "score_covariance", curvature = "observed")
  expect_identical(out$provenance$curvature, "observed")
  expect_identical(out$provenance$weight, "score_covariance")
})

test_that("weight changes the statistic, score_covariance only the reference", {
  base <- semtest(semtest_restricted, semtest_fit)
  raw <- semtest(semtest_restricted, semtest_fit, score_covariance = "uncentered")
  expected <- semtest(semtest_restricted, semtest_fit, weight = "expected_curvature")
  expect_equal(raw$statistic, base$statistic)
  expect_false(isTRUE(all.equal(raw$eigenvalues, base$eigenvalues)))
  expect_false(isTRUE(all.equal(expected$statistic, base$statistic)))
})

test_that("the score covariance is always empirical", {
  # A model-based "expected" level assumed the information equality, so it was
  # the non-robust test, whose p-value pvalues(x, "STD") already gives.
  expect_error(inference_control(score_covariance = "expected"),
               class = "semTests_error_invalid_inference")
  expect_error(semtest(semtest_fit, score_covariance = "expected"),
               class = "semTests_error_invalid_inference")
  out <- semtest(semtest_fit)
  expect_equal(as.numeric(pvalues(out, "STD")),
               stats::pchisq(out$statistic, out$df, lower.tail = FALSE))
})

test_that("both empirical score covariances work for ML and FIML", {
  d <- lavaan::HolzingerSwineford1939
  set.seed(3)
  d$x2[sample(nrow(d), 60)] <- NA
  d$x4[sample(nrow(d), 40)] <- NA
  fiml <- lavaan::cfa(semtest_model, d, missing = "ml", meanstructure = TRUE)
  fiml_restricted <- lavaan::cfa(
    paste(semtest_model, "a == b", "b == c", sep = "\n"), d,
    missing = "ml", meanstructure = TRUE
  )
  cases <- list(
    ml = list(semtest_fit), fiml = list(fiml),
    ml_nested = list(semtest_restricted, semtest_fit),
    fiml_nested = list(fiml_restricted, fiml)
  )
  for (case in names(cases)) for (covariance in c("centered", "uncentered")) {
    info <- paste(case, covariance)
    out <- do.call(semtest, c(cases[[case]], list(score_covariance = covariance)))
    expect_identical(out$provenance$score_covariance, covariance, info = info)
    expect_true(is.finite(out$p_value), info = info)
    # An empirical reference is not the flat spectrum of the non-robust test.
    expect_gt(max(abs(out$eigenvalues - 1)), 1e-3, label = info)
  }
})
