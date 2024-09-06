test_that("pvalues_one return p_values in (0, 1)", {
  pvals <- pvalues_(object, eba = 2, peba = 2, unbiased = 2, trad = "sb", pols = 2)
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})

# test_that("pvalues works with NULL", {
#   expect_error(pvalues(object, tests = NULL, trad = NULL, eba = NULL, peba = NULL, pols = NULL))
#   expect_equal(length(pvalues(object, tests = NULL, trad = NULL, eba = NULL, peba = 1, pols = NULL)), 2)
#   expect_equal(length(pvalues(object, tests = NULL, trad = NULL, eba = 1, peba = NULL, pols = 2)), 4)
#   expect_equal(length(pvalues(object, tests = NULL, trad = c("std", "sf", "ss", "sb"), eba = NULL, peba = NULL, pols = NULL)), 8)
# }, skip )
