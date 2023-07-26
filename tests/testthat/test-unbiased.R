model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

test_that("unbiased and biased not the same for object", {
  ugamma_1 <- ugamma_non_nested(object, unbiased = TRUE)
  ugamma_2 <- ugamma_non_nested(object, unbiased = FALSE)
  expect_true(sum(abs(ugamma_1 - ugamma_2)) > 1)

  ugamma_1 <- ugamma_non_nested(m1, unbiased = TRUE)
  ugamma_2 <- ugamma_non_nested(m1, unbiased = FALSE)
  expect_true(sum(abs(ugamma_1 - ugamma_2)) > 1)
})


test_that("unbiased and biased not the same for m0 m1", {
  ugamma_1 <- ugamma_nested(m0, m1, unbiased = TRUE)
  ugamma_2 <- ugamma_nested(m0, m1, unbiased = FALSE)
  expect_true(sum(abs(ugamma_1 - ugamma_2)) > 1)

  ugamma_1 <- ugamma_non_nested(m1, unbiased = TRUE)
  ugamma_2 <- ugamma_non_nested(m1, unbiased = FALSE)
  expect_true(sum(abs(ugamma_1 - ugamma_2)) > 1)
})
