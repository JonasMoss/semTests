test_that("every single-model test type returns a p-value in [0, 1]", {
  types <- c(
    "STD_RLS", "SB_RLS", "SS_RLS", "SF_RLS", "ALL_RLS", "PALL_RLS",
    "EBA2_RLS", "pEBA4_RLS", "pOLS2_RLS", "STD_ML", "SB_ML"
  )
  for (type in types) {
    p <- pvalues(object, type)
    expect_length(p, 1L)
    expect_true(all(p >= 0 & p <= 1), info = type)
  }
})

test_that("requesting several tests at once returns one p-value per test", {
  p <- pvalues(object, c("SB_RLS", "pEBA4_RLS", "pOLS2_RLS"))
  expect_length(p, 3L)
  expect_equal(names(p), c("sb_rls", "peba4_rls", "pols2_rls"))
})

test_that("pvalues_nested default (penalized, all eigenvalues) is valid", {
  p <- pvalues_nested(m0_no_groups, m1_no_groups)
  expect_length(p, 1L)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("nested p-values can be requested via direct method arguments", {
  p <- pvalues_nested_internal(m0_no_groups, m1_no_groups, tests = NULL, peba = 2)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("extras = TRUE appends chi-squares and eigenvalues", {
  base <- pvalues_internal(object, tests = NULL, peba = 4, chisq = "rls")
  extra <- pvalues_internal(object, tests = NULL, peba = 4, chisq = "rls", extras = TRUE)
  expect_gt(length(extra), length(base))
  expect_true("lambda" %in% names(extra))
})

test_that("extras = TRUE works with the unbiased gamma estimator", {
  extra <- pvalues_internal(object,
    tests = NULL, peba = 4, unbiased = 2, chisq = "rls", extras = TRUE
  )
  expect_gt(length(extra), 0L)
})
