# The nested restriction matrix (get_a_matrix). "delta" evaluates both
# Jacobians at the restricted solution and is then exact for parameter-map
# nestings, so it must agree with lavaan's exact restriction. "legacy" is the
# semTests 1.0.0 construction, which equals lavaan's default delta method.

hs <- lavaan::HolzingerSwineford1939
three <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
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
  equality = list(
    h1 = three,
    h0 = "
      visual  =~ x1 + x2 + x3
      textual =~ x4 + a*x5 + a*x6
      speed   =~ x7 + x8 + x9
    ",
    group = NULL, group.equal = NULL
  ),
  # H0 omits the path. lavaan's exact restriction needs the path present (fixed
  # at zero), so its reference fit is written that way; the model is the same.
  dropped_path = list(
    h1 = paste(three, "visual =~ x9"),
    h0 = three,
    h0_lavaan = paste(three, "visual =~ 0*x9"),
    group = NULL, group.equal = NULL
  ),
  two_group_loadings = list(
    h1 = three, h0 = three,
    group = "school", group.equal = "loadings"
  )
)
fit_pair <- function(nesting, h0 = nesting$h0) {
  m1 <- lavaan::cfa(nesting$h1, hs, estimator = "MLM", group = nesting$group)
  m0 <- lavaan::cfa(h0, hs,
    estimator = "MLM", group = nesting$group,
    group.equal = if (is.null(nesting$group.equal)) "" else nesting$group.equal
  )
  list(m0 = m0, m1 = m1)
}

test_that("delta restriction equals lavaan's exact restriction", {
  for (nm in names(nestings)) {
    fits <- fit_pair(nestings[[nm]])
    h0_ref <- if (is.null(nestings[[nm]]$h0_lavaan)) nestings[[nm]]$h0 else nestings[[nm]]$h0_lavaan
    ref <- fit_pair(nestings[[nm]], h0_ref)
    expect_equal(semtests_sb_p(fits$m0, fits$m1, "delta"),
      lavaan_sb_p(ref$m1, ref$m0, "exact"),
      tolerance = 1e-7, info = nm
    )
  }
})

test_that("legacy restriction reproduces lavaan's default delta method", {
  for (nm in names(nestings)) {
    fits <- fit_pair(nestings[[nm]])
    expect_equal(semtests_sb_p(fits$m0, fits$m1, "legacy"),
      lavaan_sb_p(fits$m1, fits$m0, "delta"),
      tolerance = 1e-7, info = nm
    )
  }
})

test_that("delta and exact restrictions agree where both apply", {
  fits <- fit_pair(nestings$equality)
  expect_equal(semtests_sb_p(fits$m0, fits$m1, "delta"),
    semtests_sb_p(fits$m0, fits$m1, "exact"),
    tolerance = 1e-8
  )
})

test_that("the unrestricted fit moves to the restricted solution", {
  fits <- fit_pair(nestings$dropped_path)
  moved <- h1_at_h0_solution(fits$m1, fits$m0)
  expect_false(is.null(moved))
  expect_equal(
    lavaan::lav_model_implied(moved@Model)$cov,
    lavaan::lav_model_implied(fits$m0@Model)$cov,
    tolerance = 1e-10
  )
})

test_that("delta falls back to legacy, with a warning, when H0 cannot be placed in H1", {
  # H1 identifies the factors by their variances, H0 by marker loadings, so the
  # parameter keys do not locate H0's solution in H1.
  m1 <- lavaan::cfa(three, hs, estimator = "MLM", std.lv = TRUE)
  m0 <- lavaan::cfa("
      visual  =~ x1 + x2 + x3
      textual =~ x4 + a*x5 + a*x6
      speed   =~ x7 + x8 + x9
    ", hs, estimator = "MLM")
  expect_null(h1_at_h0_solution(m1, m0))
  expect_warning(
    p <- semtests_sb_p(m0, m1, "delta"),
    class = "semTests_warning_legacy_restriction"
  )
  expect_equal(p, semtests_sb_p(m0, m1, "legacy"), tolerance = 1e-12)
})

test_that("FIML delta restriction equals the exact restriction", {
  data <- hs
  set.seed(21)
  data$x1[sample(nrow(data), 30L)] <- NA
  m1 <- lavaan::cfa(three, data, missing = "fiml", estimator = "MLR")
  m0 <- lavaan::cfa(nestings$equality$h0, data, missing = "fiml", estimator = "MLR")
  df <- as.integer(lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df"))
  expect_equal(
    fiml_lambdas_nested(m0, m1, df, "delta", "observed")$ug_biased,
    fiml_lambdas_nested(m0, m1, df, "exact", "observed")$ug_biased,
    tolerance = 1e-8
  )
})
