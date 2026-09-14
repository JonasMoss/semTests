test_that("tiled raw projection matches independently materialized moment covariance", {
  set.seed(331)
  X <- matrix(rexp(93), 31, 3)
  for (means in c(FALSE, TRUE)) {
    centered <- sweep(X, 2, colMeans(X), "-")
    Z <- do.call(cbind, lapply(seq_len(3), function(j) {
      do.call(cbind, lapply(j:3, function(k) centered[, j] * centered[, k]))
    }))
    Z <- sweep(Z, 2, colMeans(Z), "-")
    if (means) Z <- cbind(centered, Z)
    F <- matrix(rnorm(ncol(Z) * 2), ncol(Z), 2)
    expected <- crossprod(F, (crossprod(Z) / nrow(Z)) %*% F)
    for (tile in c(1L, 7L, 1000L)) {
      expect_equal(raw_moment_projection(X, F, means, tile), expected,
                   tolerance = 1e-12)
    }
    expect_equal(raw_moment_projection(X, F[, 1, drop = FALSE], means, 7L),
                 expected[1, 1, drop = FALSE], tolerance = 1e-12)
    expect_equal(raw_moment_projection(sweep(X, 2, c(100, -20, 5), "+"), F, means),
                 expected, tolerance = 1e-12)
  }
  expect_error(raw_moment_projection(X, matrix(1, 6, 1), tile_rows = 0), "positive")
  X[1, 1] <- NA
  expect_error(raw_moment_projection(X, matrix(1, 6, 1)), "finite complete data")
})

test_that("raw and cached projection agree across unequal groups and mean blocks", {
  for (fit in list(object, m0, m1, m0_no_groups)) {
    F <- model_effective_delta(fit)
    expect_true(can_project_raw_moments(fit))
    expect_equal(projected_gamma(fit, F, source = "raw", tile_rows = 17L),
                 projected_gamma(fit, F, source = "materialized"), tolerance = 1e-10)
  }
})

test_that("uncached single and nested ML spectra never request a full Gamma", {
  full_single <- single_model_reduced_lambdas(object, fit_df(object))
  df <- fit_df(m0_no_groups) - fit_df(m1_no_groups)
  full_nested <- lambdas_nested(m0_no_groups, m1_no_groups, df = df)
  single <- object
  unrestricted <- m1_no_groups
  single@SampleStats@NACOV <- list(NULL)
  unrestricted@SampleStats@NACOV <- list(NULL)
  local_mocked_bindings(gamma_from_lavaan = function(...) stop("full Gamma requested"))
  expect_equal(single_model_reduced_lambdas(single, fit_df(single)), full_single,
               tolerance = 1e-10)
  expect_equal(lambdas_nested(m0_no_groups, unrestricted, df = df), full_nested,
               tolerance = 1e-10)
})

test_that("ML fits without precomputed NACOV work directly from observations", {
  fit <- lavaan::cfa(hs_model, lavaan::HolzingerSwineford1939,
                     estimator = "ML", se = "none", test = "standard")
  expect_true(all(vapply(fit@SampleStats@NACOV, is.null, logical(1))))
  expected <- lavaan_lambdas(fit, fit_df(fit))
  local_mocked_bindings(gamma_from_lavaan = function(...) stop("full Gamma requested"))
  expect_equal(single_model_reduced_lambdas(fit, fit_df(fit)), expected,
               tolerance = 1e-10)
  expect_true(is.finite(as.numeric(pvalues(fit, "SB_ML"))))
})

test_that("cached or explicitly supplied Gamma is preserved", {
  F <- model_effective_delta(object)
  expected <- projected_gamma(object, F, source = "materialized")
  custom <- lavaan::cfa(hs_model, lavaan::HolzingerSwineford1939,
                        estimator = "MLM", NACOV = 2 * gamma_from_lavaan(object)[[1]])
  expect_false(can_project_raw_moments(custom))
  local_mocked_bindings(raw_moment_projection = function(...) stop("raw rows requested"))
  expect_equal(projected_gamma(object, F), expected)
  expect_equal(projected_gamma(custom, F), 2 * expected, tolerance = 1e-10)
  expect_error(projected_gamma(custom, F, source = "raw"), "not supported")
})

test_that("unverified first-stage conventions cannot enter the raw route", {
  for (option in c("gamma.n.minus.one", "gamma.unbiased")) {
    fit <- object
    fit@Options[[option]] <- TRUE
    expect_false(can_project_raw_moments(fit))
  }
  fit <- object
  fit@Options$likelihood <- "wishart"
  expect_false(can_project_raw_moments(fit))
  fit <- object
  fit@Data@weights <- list(rep(1, fit@SampleStats@nobs[[1]]))
  expect_false(can_project_raw_moments(fit))
  expect_false(can_project_raw_moments(m1_))
})
