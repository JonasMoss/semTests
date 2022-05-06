test_that("pvalues_two return p_values in (0, 1)", {
  pvals <- pvalues_two(m0, m1)
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})

# test_that("pvalues_one return p_values in (0, 1)", {
#   pvals <- pvalues_one(object)
#   expect_true(all(pvals < 1))
#   expect_true(all(pvals > 0))
# })
