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
