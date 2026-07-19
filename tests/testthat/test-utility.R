tests <- c("SB_UG_RLS", "pEBA2_UG_RLS", "EBA4_RLS", "pEBA6_RLS", "pOLS2_UG_ML")
options <- lapply(tests, \(test) split_input(test))
result <- unlist(lapply(options, \(option) {
  do.call(compute_pvalues, c(list(m0 = object), option))
}))

test_that("split_input applied to p_values gives correct names", {
  expect_equal(names(result), tolower(tests))
})

test_that("pvalues dispatches parsed test specifications to the numeric engine", {
  actual <- pvalues(object, tests)
  expect_equal(as.numeric(actual), as.numeric(result))
  expect_equal(names(actual), names(result))
})
