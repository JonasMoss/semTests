# The production nested path computes the reference spectrum from the reduced
# m x m matrix (lambdas_nested). These tests pin it to the full q x q UGamma
# eigendecomposition (ugamma_nested_reference + base eigen), which must agree.

ndf <- function(m0, m1) {
  as.numeric(lavaan::fitMeasures(m0, "df")) - as.numeric(lavaan::fitMeasures(m1, "df"))
}

# Full reference: top-df eigenvalues of the full q x q UGamma, per gamma estimate.
full_spectrum <- function(m0, m1, unbiased, df) {
  ug_list <- ugamma_nested_reference(m0, m1, unbiased)
  lapply(ug_list, function(ug) {
    sort(Re(eigen(ug, only.values = TRUE)$values), decreasing = TRUE)[seq_len(df)]
  })
}

expect_spectrum_matches <- function(m0, m1, unbiased) {
  df <- ndf(m0, m1)
  red <- lambdas_nested(m0, m1, unbiased, df)
  full <- full_spectrum(m0, m1, unbiased, df)
  expect_equal(length(red), length(full))
  for (j in seq_along(red)) {
    expect_equal(sort(red[[j]]), sort(full[[j]]), tolerance = 1e-8)
  }
}

test_that("reduced 2000 spectrum matches full UGamma (no groups, biased + unbiased)", {
  skip_on_cran() # cross-route numeric-equivalence; not portable across BLAS/LAPACK
  expect_spectrum_matches(m0_no_groups, m1_no_groups, 1)
  expect_spectrum_matches(m0_no_groups, m1_no_groups, 2)
})

test_that("reduced 2000 spectrum matches full UGamma (multi-group metric invariance)", {
  skip_on_cran() # cross-route numeric-equivalence; not portable across BLAS/LAPACK
  expect_spectrum_matches(m0, m1, 1)
  expect_spectrum_matches(m0, m1, 2)
})

test_that("reduced spectrum matches full UGamma under a general (==) equality constraint", {
  skip_on_cran() # cross-route numeric-equivalence; not portable across BLAS/LAPACK
  # `a*` labels create *simple* equalities (ceq.simple); an explicit `==` makes a
  # general equality constraint (eq.constraints), which routes through the
  # eq.constraints.K mapping in both the reduced and the full reference paths.
  HS <- lavaan::HolzingerSwineford1939
  m1 <- lavaan::cfa("visual =~ x1 + x2 + x3
                     textual =~ x4 + x5 + x6
                     speed   =~ x7 + x8 + x9", HS, estimator = "MLM")
  m0 <- lavaan::cfa("visual =~ x1 + v2*x2 + v3*x3
                     textual =~ x4 + x5 + x6
                     speed   =~ x7 + x8 + x9
                     v2 == v3", HS, estimator = "MLM")
  expect_true(m0@Model@eq.constraints) # general constraint, not ceq.simple
  expect_spectrum_matches(m0, m1, 1)
  expect_spectrum_matches(m0, m1, 2)
})

test_that("nested continuous tests support an already-constrained H1", {
  skip_on_cran() # cross-route numeric-equivalence; not portable across BLAS/LAPACK
  HS <- lavaan::HolzingerSwineford1939
  h1 <- "
    visual  =~ x1 + a*x2 + b*x3
    textual =~ x4 + c*x5 + d*x6
    speed   =~ x7 + x8 + x9
    a == b
  "
  h0 <- paste(h1, "c == d", sep = "\n")
  m1 <- lavaan::cfa(h1, HS, estimator = "MLM")
  m0 <- lavaan::cfa(h0, HS, estimator = "MLM")

  expect_true(m1@Model@eq.constraints)
  expect_spectrum_matches(m0, m1, 1)
  expect_spectrum_matches(m0, m1, 2)
  expect_true(all(pvalues_nested(m0, m1) >= 0))
})

test_that("nested continuous tests support equality mixed with inequality", {
  HS <- lavaan::HolzingerSwineford1939
  h1 <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  "
  h0 <- "
    visual  =~ x1 + a*x2 + b*x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
    a == b
    a > .1
  "
  m1 <- lavaan::cfa(h1, HS, estimator = "MLM")
  m0 <- lavaan::cfa(h0, HS, estimator = "MLM")

  expect_false(m0@Model@eq.constraints)
  expect_gt(nrow(m0@Model@ceq.JAC), 0L)
  A <- get_a_matrix(m1, m0)
  expect_equal(nrow(A), ndf(m0, m1))
  expect_true(all(is.finite(pvalues_nested(m0, m1))))
})

test_that("restriction maps retain lavaan simple-constraint compatibility", {
  legacy <- m1_no_groups
  npar <- legacy@Model@nx.free
  legacy@Model@eq.constraints <- FALSE
  legacy@Model@ceq.simple.only <- TRUE
  legacy@Model@ceq.simple.K <- diag(npar)[, -npar, drop = FALSE]

  A <- get_a_matrix(m1_no_groups, legacy)
  expect_equal(dim(A), c(1L, npar))
})

test_that("reduced spectrum matches full UGamma on a big model (q = 210)", {
  skip_on_cran()
  set.seed(1)
  p <- 20
  ind <- paste0("y", 1:p)
  pop <- paste("f =~", paste(sprintf("0.7*%s", ind), collapse = " + "))
  dat <- lavaan::simulateData(pop, sample.nobs = 600)
  m1_big <- lavaan::cfa(paste("f =~", paste(ind, collapse = " + ")),
    dat,
    estimator = "MLM"
  )
  tied <- c(ind[1], sprintf("a*%s", ind[2:10]), ind[11:p]) # 8 loading equalities
  m0_big <- lavaan::cfa(paste("f =~", paste(tied, collapse = " + ")),
    dat,
    estimator = "MLM"
  )
  q <- nrow(gamma_matrices(m1_big, 1, m0_big)[["ug_biased"]])
  expect_gt(q, 200L) # genuinely big: q = 210, while m = 8
  expect_spectrum_matches(m0_big, m1_big, 1)
  expect_spectrum_matches(m0_big, m1_big, 2)

  # the reduced path must not be slower than forming the full q x q UGamma
  t_red <- system.time(
    lambdas_nested(m0_big, m1_big, 1, ndf(m0_big, m1_big))
  )[["elapsed"]]
  t_full <- system.time(
    full_spectrum(m0_big, m1_big, 1, ndf(m0_big, m1_big))
  )[["elapsed"]]
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
    dat,
    estimator = "MLM"
  )
  tied <- c(ind[1], sprintf("a*%s", ind[2:10]), ind[11:p])
  m0_big <- lavaan::cfa(paste("f =~", paste(tied, collapse = " + ")),
    dat,
    estimator = "MLM"
  )
  pv <- pvalues_nested(m0_big, m1_big)
  expect_true(all(pv >= 0 & pv <= 1))
})
