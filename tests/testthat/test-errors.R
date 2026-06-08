test_that("pvalues errors when no p-values are requested", {
  expect_error(pvalues(object, tests = NULL), "provide some p-values")
})

test_that("pvalues_nested errors when no p-values are requested", {
  expect_error(
    pvalues_nested(m0_no_groups, m1_no_groups, tests = NULL),
    "provide some p-values"
  )
})

test_that("pvalues_nested errors when the models have equal degrees of freedom", {
  expect_error(
    pvalues_nested(m1_no_groups, m1_no_groups),
    "same degree of freedom"
  )
})

test_that("pvalues_nested now supports continuous non-ML estimators (GLS)", {
  p <- pvalues_nested(m0_, m1_)            # GLS fits from setup.R
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
})
