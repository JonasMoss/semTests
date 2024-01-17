model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

test_that("unbiased and biased not the same for object", {
  ugamma_1 <- ugamma_no_groups(object, unbiased = 1)$ug_biased
  ugamma_2 <- ugamma_no_groups(object, unbiased = 2)$ug_unbiased
  expect_true(sum(abs(ugamma_1 - ugamma_2)) > 1)
})
