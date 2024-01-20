tests <- c("SB_UG_RLS", "pEBA2_UG_RLS", "EBA4_RLS", "pEBA6_RLS", "pOLS2_UG")
options <- sapply(tests, \(test) split_input(test))

test_that("split_input looks reasonable", {
  expect_equal(colnames(options), tests)
  expect_equal(dim(options), c(6, 5))
})

options <- lapply(tests, \(test) split_input(test))
result <- sapply(options, \(option) do.call(pvalues_one, c(object,option)))
test_that("split_input applied to p_values gives correct names", {
  expect_equal(names(result), tolower(tests))
})

test_that("using 'tests' in pvalues yield the correct results.", {
  expect_equal(pvalues(object, tests[1]),
               pvalues(object,
                       tests = NULL,
                       trad = "sb",
                       peba = NULL,
                       eba = NULL,
                       pols = NULL,
                       unbiased = 2,
                       chisq = "rls"))


  expect_equal(pvalues(object, tests[2]),
               pvalues(object,
                       tests = NULL,
                       trad = NULL,
                       peba = 2,
                       eba = NULL,
                       pols = NULL,
                       unbiased = 2,
                       chisq = "rls"))


  expect_equal(pvalues(object, tests[3]),
               pvalues(object,
                       tests = NULL,
                       trad = NULL,
                       peba = NULL,
                       eba = 4,
                       pols = NULL,
                       unbiased = 1,
                       chisq = "rls"))


  expect_equal(pvalues(object, tests[4]),
               pvalues(object,
                       tests = NULL,
                       trad = NULL,
                       peba = 6,
                       eba = NULL,
                       pols = NULL,
                       unbiased = 1,
                       chisq = "rls"))

  expect_equal(pvalues(object, tests[5]),
               pvalues(object,
                       tests = NULL,
                       trad = NULL,
                       peba = NULL,
                       eba = NULL,
                       pols = 2,
                       unbiased = 2,
                       chisq = "ml"))
})
