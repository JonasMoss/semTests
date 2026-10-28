# The nested restriction matrix (get_a_matrix). Both constructions carry
# lavaan's names: "delta" (the default) is lavaan's default moment-tangent map,
# and "exact" is lavaan's literal parameter restriction map.

hs <- lavaan::HolzingerSwineford1939
three <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
equal_loadings <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + a*x5 + a*x6
  speed   =~ x7 + x8 + x9
"

lavaan_sb_p <- function(m1, m0, A.method) {
  lavaan::lavTestLRT(m1, m0,
    method = "satorra.2000", test = "satorra.bentler",
    scaled_shifted = FALSE, A.method = A.method
  )[["Pr(>Chisq)"]][2]
}
semtests_sb_p <- function(m0, m1, A.method) {
  as.numeric(unlist(pvalues_nested(m0, m1, tests = "SB_ML", A.method = A.method)))
}

nestings <- list(
  equality = list(h1 = three, h0 = equal_loadings, group = NULL, group.equal = NULL),
  dropped_path = list(
    h1 = paste(three, "visual =~ x9"), h0 = three,
    group = NULL, group.equal = NULL
  ),
  two_group_loadings = list(
    h1 = three, h0 = three,
    group = "school", group.equal = "loadings"
  )
)
fit_pair <- function(nesting, ...) {
  m1 <- lavaan::cfa(nesting$h1, hs, estimator = "MLM", group = nesting$group, ...)
  m0 <- lavaan::cfa(nesting$h0, hs,
    estimator = "MLM", group = nesting$group,
    group.equal = if (is.null(nesting$group.equal)) "" else nesting$group.equal,
    ...
  )
  list(m0 = m0, m1 = m1)
}

test_that("delta restriction reproduces lavaan's default delta method", {
  skip_on_cran() # cross-implementation agreement
  for (nm in names(nestings)) {
    fits <- fit_pair(nestings[[nm]])
    expect_equal(semtests_sb_p(fits$m0, fits$m1, "delta"),
      lavaan_sb_p(fits$m1, fits$m0, "delta"),
      tolerance = 1e-7, info = nm
    )
  }
})

test_that("exact restriction reproduces lavaan's exact method for complete data", {
  skip_on_cran() # cross-implementation agreement
  for (nm in c("equality", "two_group_loadings")) {
    fits <- fit_pair(nestings[[nm]])
    expect_equal(semtests_sb_p(fits$m0, fits$m1, "exact"),
      lavaan_sb_p(fits$m1, fits$m0, "exact"),
      tolerance = 1e-7, info = nm
    )
  }
})

test_that("exact restriction refuses models with different parameter sets", {
  fits <- fit_pair(nestings$dropped_path)
  expect_error(
    semtests_sb_p(fits$m0, fits$m1, "exact"),
    "same underlying full free-parameter set"
  )
})

test_that("delta restriction does not depend on the factor identification", {
  # Marker loadings and std.lv = TRUE are the same models on the same data, so
  # only optimizer noise may separate them (about 5e-6 relative here). The
  # exact restriction, which depends on the identification, moves this p-value
  # from .0114 to .0101.
  marker <- fit_pair(nestings$equality)
  standardized <- fit_pair(nestings$equality, std.lv = TRUE)
  expect_equal(semtests_sb_p(marker$m0, marker$m1, "delta"),
    semtests_sb_p(standardized$m0, standardized$m1, "delta"),
    tolerance = 1e-4
  )
})

test_that("A.method is recorded for complete-data nested tests", {
  fits <- fit_pair(nestings$equality)
  for (A.method in c("delta", "exact")) {
    result <- pvalues_nested(fits$m0, fits$m1, tests = "SB_ML", A.method = A.method)
    expect_identical(attr(result, "semtests")$A.method, A.method)
  }
})
