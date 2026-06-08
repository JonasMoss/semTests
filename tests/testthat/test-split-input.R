test_that("split_input parses the two-part unbiased form", {
  opt <- split_input("sb_ug")
  expect_equal(opt$unbiased, 2)
  expect_equal(opt$trad, "sb")
  expect_equal(opt$chisq, "rls")
})

test_that("split_input parses the two-part chi-square form", {
  expect_equal(split_input("sb_ml")$chisq, "ml")
  expect_equal(split_input("sb_rls")$chisq, "rls")
})

test_that("split_input rejects unrecognised test strings", {
  expect_error(split_input("garbage"), "Invalid")
})

test_that("default() returns 2 for empty input and the number otherwise", {
  expect_equal(default(""), 2)
  expect_equal(default("5"), 5)
})

test_that("nanull() maps NA to NULL and passes other values through", {
  expect_null(nanull(NA))
  expect_equal(nanull(3), 3)
})

test_that("pols_pvalue handles a single eigenvalue", {
  expect_true(is.finite(suppressWarnings(pols_pvalue(22, 1, 2))))
})
