# The production nested path computes the reference spectrum from the reduced
# m x m matrix (lambdas_nested). These tests pin it to the full q x q UGamma
# eigendecomposition (ugamma_nested + base eigen), which must agree exactly.

ndf <- function(m0, m1) {
  as.numeric(lavaan::fitMeasures(m0, "df")) - as.numeric(lavaan::fitMeasures(m1, "df"))
}

# Full reference: top-df eigenvalues of the full q x q UGamma, per gamma estimate.
full_spectrum <- function(m0, m1, method, unbiased, df) {
  ug_list <- ugamma_nested(m0, m1, method, unbiased)
  lapply(ug_list, function(ug) {
    sort(Re(eigen(ug, only.values = TRUE)$values), decreasing = TRUE)[seq_len(df)]
  })
}

expect_spectrum_matches <- function(m0, m1, method, unbiased) {
  df <- ndf(m0, m1)
  red <- lambdas_nested(m0, m1, method, unbiased, df)
  full <- full_spectrum(m0, m1, method, unbiased, df)
  expect_equal(length(red), length(full))
  for (j in seq_along(red)) {
    expect_equal(sort(red[[j]]), sort(full[[j]]), tolerance = 1e-8)
  }
}

test_that("reduced 2000 spectrum matches full UGamma (no groups, biased + unbiased)", {
  expect_spectrum_matches(m0_no_groups, m1_no_groups, "2000", 1)
  expect_spectrum_matches(m0_no_groups, m1_no_groups, "2000", 2)
})

test_that("reduced 2000 spectrum matches full UGamma (multi-group metric invariance)", {
  expect_spectrum_matches(m0, m1, "2000", 1)
  expect_spectrum_matches(m0, m1, "2000", 2)
})

test_that("2001 spectrum matches full UGamma top-df (no reduction, base eigen)", {
  expect_spectrum_matches(m0_no_groups, m1_no_groups, "2001", 1)
  expect_spectrum_matches(m0_no_groups, m1_no_groups, "2001", 2)
})

test_that("reduced spectrum matches full UGamma on a big model (q = 210)", {
  skip_on_cran()
  set.seed(1)
  p <- 20
  ind <- paste0("y", 1:p)
  pop <- paste("f =~", paste(sprintf("0.7*%s", ind), collapse = " + "))
  dat <- lavaan::simulateData(pop, sample.nobs = 600)
  m1_big <- lavaan::cfa(paste("f =~", paste(ind, collapse = " + ")),
                        dat, estimator = "MLM")
  tied <- c(ind[1], sprintf("a*%s", ind[2:10]), ind[11:p]) # 8 loading equalities
  m0_big <- lavaan::cfa(paste("f =~", paste(tied, collapse = " + ")),
                        dat, estimator = "MLM")
  q <- nrow(gamma(m1_big, 1, m0_big)[["ug_biased"]])
  expect_gt(q, 200L) # genuinely big: q = 210, while m = 8
  expect_spectrum_matches(m0_big, m1_big, "2000", 1)
  expect_spectrum_matches(m0_big, m1_big, "2000", 2)

  # the reduced path must not be slower than forming the full q x q UGamma
  t_red <- system.time(lambdas_nested(m0_big, m1_big, "2000", 1, ndf(m0_big, m1_big)))[["elapsed"]]
  t_full <- system.time(full_spectrum(m0_big, m1_big, "2000", 1, ndf(m0_big, m1_big)))[["elapsed"]]
  expect_lte(t_red, t_full + 0.05)
})

test_that("end-to-end pvalues_nested runs on the big model", {
  skip_on_cran()
  set.seed(2)
  p <- 20
  ind <- paste0("y", 1:p)
  pop <- paste("f =~", paste(sprintf("0.7*%s", ind), collapse = " + "))
  dat <- lavaan::simulateData(pop, sample.nobs = 600)
  m1_big <- lavaan::cfa(paste("f =~", paste(ind, collapse = " + ")),
                        dat, estimator = "MLM")
  tied <- c(ind[1], sprintf("a*%s", ind[2:10]), ind[11:p])
  m0_big <- lavaan::cfa(paste("f =~", paste(tied, collapse = " + ")),
                        dat, estimator = "MLM")
  pv <- pvalues_nested(m0_big, m1_big)
  expect_true(all(pv >= 0 & pv <= 1))
})
