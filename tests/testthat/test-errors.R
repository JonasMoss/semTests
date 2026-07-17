test_that("pvalues errors when no p-values are requested", {
  expect_error(
    pvalues(object, tests = NULL),
    "tests.*at least one",
    class = "semTests_error_invalid_tests"
  )
})

test_that("pvalues_nested errors when no p-values are requested", {
  expect_error(
    pvalues_nested(m0_no_groups, m1_no_groups, tests = NULL),
    "tests.*at least one",
    class = "semTests_error_invalid_tests"
  )
})

test_that("pvalues_nested errors when the models have equal degrees of freedom", {
  expect_error(
    pvalues_nested(m1_no_groups, m1_no_groups),
    "different degrees of freedom",
    class = "semTests_error_incompatible_models"
  )
})

test_that("pvalues_nested now supports continuous non-ML estimators (GLS)", {
  p <- pvalues_nested(m0_, m1_)            # GLS fits from setup.R
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
})

test_that("public test specifications fail early and clearly", {
  invalid <- list(
    character(),
    NA_character_,
    "",
    1,
    "PEBA0",
    "PEBA2.5",
    "PEBA1e20",
    "POLS0",
    "POLS-1",
    "SB_NOPE",
    "SB_UG_NOPE",
    c("SB", NA_character_)
  )
  for (tests in invalid) {
    expect_error(
      pvalues(object, tests = tests),
      class = "semTests_error_invalid_tests",
      info = paste(deparse(tests), collapse = "")
    )
  }
  expect_error(
    pvalues(object, tests = "PEBA1000"),
    "more blocks than the test degrees of freedom",
    class = "semTests_error_invalid_tests"
  )
})

test_that("fit-quality failures have stable public conditions", {
  nonconverged <- object
  nonconverged@optim$converged <- FALSE
  expect_error(
    pvalues(nonconverged, "SB_ML"),
    "did not converge",
    class = "semTests_error_nonconverged"
  )

  inadmissible <- object
  variance_row <- which(
    inadmissible@ParTable$op == "~~" &
      inadmissible@ParTable$lhs == inadmissible@ParTable$rhs &
      inadmissible@ParTable$lhs == "x1"
  )
  inadmissible@ParTable$est[variance_row[1L]] <- -1
  expect_warning(
    check_fit_quality(inadmissible),
    "inadmissible",
    class = "semTests_warning_inadmissible"
  )
})

test_that("nested fits must be comparable before a difference test is run", {
  different_data <- lavaan::HolzingerSwineford1939
  different_data$x1 <- different_data$x1 +
    seq_len(nrow(different_data)) / 1000
  other_data_fit <- lavaan::cfa(
    hs_model, different_data, estimator = "MLM"
  )
  expect_error(
    pvalues_nested(m0_no_groups, other_data_fit),
    "same raw data",
    class = "semTests_error_incompatible_models"
  )

  other_estimator <- lavaan::cfa(
    hs_model, lavaan::HolzingerSwineford1939, estimator = "MLR"
  )
  expect_error(
    pvalues_nested(m0_no_groups, other_estimator),
    "same requested lavaan estimator",
    class = "semTests_error_incompatible_models"
  )

  other_information <- m1_no_groups
  other_information@Options$information <- c("observed", "observed")
  expect_error(
    pvalues_nested(m0_no_groups, other_information),
    "information",
    class = "semTests_error_incompatible_models"
  )
})

test_that("reversed nested inputs are swapped with an explicit warning", {
  forward <- pvalues_nested(m0_no_groups, m1_no_groups)
  reversed <- NULL
  expect_warning(
    reversed <- pvalues_nested(m1_no_groups, m0_no_groups),
    "swapping the inputs",
    class = "semTests_warning_model_order"
  )
  expect_equal(as.numeric(forward), as.numeric(reversed))
})

test_that("method 2001 never silently falls back to method 2000", {
  testthat::local_mocked_bindings(
    lambdas_nested = function(...) list(ug_biased = c(1, -0.1)),
    .package = "semTests"
  )
  expect_error(
    pvalues_nested(
      m0_no_groups, m1_no_groups,
      method = "2001", tests = "SB_ML"
    ),
    "retry with `method = \"2000\"`",
    class = "semTests_error_unstable_spectrum"
  )
})
