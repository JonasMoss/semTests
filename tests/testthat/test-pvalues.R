test_that("pvalues_two return p_values in (0, 1)", {
  pvals <- pvalues_two(m0, m1)
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})

test_that("estimators other than ml do not work.", {
  expect_error(pvalues_two(m0_, m1_))
  expect_error(pvalues_one(m0_))
})

test_that("pvalues_one return p_values in (0, 1)", {
  pvals <- pvalues_one(m1)
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})

test_that("pvalues_two return p_values in (0, 1)", {
  expect_equal(pvalues(m1), pvalues_one(m1))
  expect_equal(pvalues_two(m0, m1), pvalues(m0, m1))
})

test_that("models with groups and equality not supported.", {
  expect_error(pvalues_one(m0))
})
