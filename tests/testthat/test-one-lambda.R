chisq <- 22
lambdas <- 1

test_that("a single lambda is handled gracefully (no error), returning a finite p-value", {
  expect_true(is.finite(peba_pvalue(chisq, lambdas, 2)))
  expect_true(is.finite(eba_pvalue(chisq, lambdas, 2)))
  expect_true(is.finite(pols_pvalue(chisq, lambdas, 2)))
})
