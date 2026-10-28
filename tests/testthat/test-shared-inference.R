test_that("batched calibrations preserve ordering, duplicates and gamma choices", {
  specs <- c("SB_ML", "STD_ML", "SS_UG_ML", "SB_ML", "PEBA2_RLS", "SB")
  for (nested in c(FALSE, TRUE)) {
    run <- if (nested) {
      function(tests) pvalues_nested(m0_no_groups, m1_no_groups, tests = tests)
    } else {
      function(tests) pvalues(object, tests)
    }
    separate <- lapply(specs, run)
    batch <- run(specs)
    expect_equal(as.numeric(batch), unlist(lapply(separate, as.numeric)))
    expect_identical(names(batch), unlist(lapply(separate, names)))
    expect_identical(attr(batch, "semtests")$requested_tests, specs)
  }
})

test_that("one call shares spectra across statistics and resolved aliases", {
  original_spectra <- compute_spectra
  original_statistic <- make_chisqs
  spectra <- statistics <- 0L
  local_mocked_bindings(
    compute_spectra = function(...) {
      spectra <<- spectra + 1L
      original_spectra(...)
    },
    make_chisqs = function(...) {
      statistics <<- statistics + 1L
      original_statistic(...)
    }
  )
  specs <- c("SB", "SS_RLS", "SB_ML", "SS_UG_ML", "SB_UG_RLS")
  pvalues(object, specs)
  expect_equal(spectra, 2L)
  expect_equal(statistics, 2L)
  # A new call, with different models, must never reuse the previous context.
  pvalues_nested(m0_no_groups, m1_no_groups, tests = specs)
  expect_equal(spectra, 4L)
  expect_equal(statistics, 4L)
})

test_that("standard-only requests avoid all spectral work", {
  local_mocked_bindings(compute_spectra = function(...) stop("spectrum requested"))
  actual <- pvalues(object, "STD_ML")
  expected <- 1 - pchisq(make_chisqs("ml", object), fit_df(object))
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-8)
  actual <- pvalues_nested(m0_no_groups, m1_no_groups, tests = "STD_ML")
  expected <- 1 - pchisq(make_chisqs("ml", m0_no_groups, m1_no_groups),
                    fit_df(m0_no_groups) - fit_df(m1_no_groups))
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-8)
  expect_error(pvalues(m1_, "STD_UG_ML"), "unbiased")
  expect_error(pvalues(m1_, "STD_RLS"), "RLS statistic")
})
