chisq <- 22
lambdas <- 1
gamma <- 2

test_that("warning, not error, on one lambda", {
  expect_warning(peba_pvalue(chisq, lambdas, 2))
  expect_warning(eba_pvalue(chisq, lambdas, 2))
  expect_warning(peba_pvalue(chisq, lambdas, 2))
})
