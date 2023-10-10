test_that("pvalues_two return p_values in (0, 1)", {
  pvals <- pvalues_two(m0, m1, eba = 2, trad = "psb")
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})

test_that("estimators other than ml do not work.", {
  expect_error(pvalues_two(m0_, m1_, eba = 2, trad = "psb"))
  # expect_error(pvalues_one(m0_, eba = 2, trad = "psb", unbiased = 2))
})

test_that("pvalues_one return p_values in (0, 1)", {
  pvals <- pvalues_one(object, eba = 2, eba_half = 2, unbiased = 2, trad = "psb")
  expect_true(all(pvals < 1))
  expect_true(all(pvals > 0))
})

test_that("pvalues_two and pvalues / pvalues_one and pvalues agree", {
  expect_equal(
    pvalues(object, eba = 2, trad = "psb", eba_half = 2, unbiased = 2),
    pvalues_one(object, eba = 2, trad = "psb",  eba_half = 2, unbiased = 2)
  )
  expect_equal(
    pvalues_two(m0, m1, eba = 2, trad = "psb", unbiased = FALSE),
    pvalues(m0, m1, eba = 2, trad = "psb", unbiased = 1)
  )
})

test_that("models with groups and equality not supported.", {
  expect_error(pvalues_one(m0, eba = c(1:2), trad = "psb", unbiased = TRUE))
})
