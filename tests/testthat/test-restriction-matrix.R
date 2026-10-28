# The nested restriction matrix (get_a_matrix). "weighted" (the default) is the
# moment-tangent map projected in the fit's weight metric. "delta" is the same
# map in the Euclidean metric, lavaan's default, and "exact" is lavaan's
# literal parameter restriction map.

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

test_that("weighted and delta restrictions do not depend on the factor identification", {
  # Marker loadings and std.lv = TRUE are the same models on the same data, so
  # only optimizer noise may separate them (about 5e-6 relative here). The
  # exact restriction, which depends on the identification, moves this p-value
  # from .0114 to .0101.
  marker <- fit_pair(nestings$equality)
  standardized <- fit_pair(nestings$equality, std.lv = TRUE)
  for (A.method in c("weighted", "delta")) {
    expect_equal(semtests_sb_p(marker$m0, marker$m1, A.method),
      semtests_sb_p(standardized$m0, standardized$m1, A.method),
      tolerance = 1e-4, info = A.method
    )
  }
})

test_that("the weighted restriction does not depend on the units of the variables", {
  # Multiplying x4 by 10 changes the delta p-value from .0072 to .0085, and
  # the weighted one (.0082) only through optimizer noise. Larger factors push
  # lavaan's own optimizer to slightly different solutions.
  rescaled <- hs
  rescaled$x4 <- 10 * rescaled$x4
  fit_scaled <- function(data, ...) {
    list(m0 = lavaan::cfa(equal_loadings, data, ...), m1 = lavaan::cfa(three, data, ...))
  }
  # lavaan fits a label shared across groups imprecisely after rescaling, so
  # the grouped pair constrains loadings with group.equal instead.
  fit_groups <- function(data) {
    list(m0 = lavaan::cfa(three, data, group = "school", group.equal = "loadings"),
         m1 = lavaan::cfa(three, data, group = "school"))
  }
  set.seed(3)
  missing <- hs
  missing$x5[sample(nrow(hs), 50)] <- NA
  missing_scaled <- missing
  missing_scaled$x4 <- 10 * missing_scaled$x4
  pairs <- list(
    complete = list(fit_scaled(hs), fit_scaled(rescaled)),
    fiml = list(fit_scaled(missing, missing = "fiml"),
                fit_scaled(missing_scaled, missing = "fiml")),
    groups = list(fit_groups(hs), fit_groups(rescaled))
  )
  for (nm in names(pairs)) {
    original <- pairs[[nm]][[1]]
    scaled <- pairs[[nm]][[2]]
    expect_equal(semtests_sb_p(scaled$m0, scaled$m1, "weighted"),
      semtests_sb_p(original$m0, original$m1, "weighted"),
      tolerance = 1e-4, info = nm
    )
  }
  original <- pairs$complete[[1]]
  scaled <- pairs$complete[[2]]
  expect_gt(abs(semtests_sb_p(scaled$m0, scaled$m1, "delta") -
                  semtests_sb_p(original$m0, original$m1, "delta")), 1e-3)
})

test_that("the weighted restriction keeps its rank under extreme units", {
  # Variances 1e8 apart make Delta' W Delta numerically singular in the raw
  # parameters; the standardized coordinates keep the single restriction.
  extreme <- hs
  extreme$x1 <- 1e4 * extreme$x1
  extreme$x4 <- extreme$x4 / 100
  # lavaan warns about the variance ratios, which is the point here.
  m1 <- suppressWarnings(lavaan::cfa(three, extreme))
  m0 <- suppressWarnings(lavaan::cfa(equal_loadings, extreme))
  expect_identical(nrow(get_a_matrix(m1, m0, "weighted")), 1L)
  expect_equal(semtests_sb_p(m0, m1, "weighted"),
    semtests_sb_p(lavaan::cfa(equal_loadings, hs), lavaan::cfa(three, hs), "weighted"),
    tolerance = 1e-4
  )
})

test_that("the weighted restriction equals delta when the weight is the identity", {
  m1 <- lavaan::cfa(three, hs, estimator = "ULS")
  m0 <- lavaan::cfa(equal_loadings, hs, estimator = "ULS")
  expect_equal(semtests_sb_p(m0, m1, "weighted"), semtests_sb_p(m0, m1, "delta"),
               tolerance = 1e-10)
})

test_that("scale helpers equilibrate information-type matrices", {
  expect_equal(positive_scale(c(4, 0, NA, -1, Inf)), c(2, 1, 1, 1, 1))
  x <- diag(c(1e8, 1, 1e-8))
  x[1, 2] <- x[2, 1] <- 1e3
  expect_equal(equilibrated_inverse(x) %*% x, diag(3), tolerance = 1e-8)
})

test_that("A.method is recorded for complete-data nested tests", {
  fits <- fit_pair(nestings$equality)
  for (A.method in c("weighted", "delta", "exact")) {
    result <- pvalues_nested(fits$m0, fits$m1, tests = "SB_ML", A.method = A.method)
    expect_identical(attr(result, "semtests")$A.method, A.method)
  }
  default <- pvalues_nested(fits$m0, fits$m1, tests = "SB_ML")
  expect_identical(attr(default, "semtests")$A.method, "weighted")
})
