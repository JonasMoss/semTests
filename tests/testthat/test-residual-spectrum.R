test_that("residual factor reproduces the full projection and ignores parameter units", {
  set.seed(9102)
  X <- matrix(rnorm(64), 8)
  H <- crossprod(X) + diag(8)
  D <- matrix(rnorm(24), 8)
  B <- residual_metric_factor(list(H), D, 5L)
  HD <- H %*% D
  U <- H - HD %*% solve(crossprod(D, HD), t(HD))
  expect_equal(tcrossprod(B), U, tolerance = 1e-8)
  scaled <- residual_metric_factor(list(H), sweep(D, 2, c(1e-8, 1, 1e8), "*"), 5L)
  expect_equal(tcrossprod(scaled), U, tolerance = 1e-8)
  expect_error(residual_metric_factor(list(H), cbind(D, D[, 1]), 4L), "rank")
  expect_null(residual_metric_factor(list(-H), D, 5L))
  expect_error(positive_symmetric_eigenvalues(diag(c(1, -1)), 2, "test"),
               "negative")
  expect_error(positive_symmetric_eigenvalues(diag(c(1, 0)), 2, "test"), "rank")
})

test_that("symmetric restriction solves agree with the general product", {
  C <- matrix(c(2, .3, .3, 1), 2)
  S <- matrix(c(3, .4, .4, 2), 2)
  expect_equal(restriction_spectrum(C, S, 2, "test"),
               ugamma_eigenvalues(solve(C, S), 2), tolerance = 1e-8)
  # Changing restriction coordinates leaves the generalized spectrum unchanged.
  T <- matrix(c(2, 1, 0, 3), 2)
  expect_equal(restriction_spectrum(crossprod(T, C %*% T),
                                    crossprod(T, S %*% T), 2, "test"),
               restriction_spectrum(C, S, 2, "test"), tolerance = 1e-8)
})

test_that("expected-information spectra match full lavaan UGamma", {
  general <- lavaan::cfa(paste(hs_model, "x2 ~~ a*x2\nx3 ~~ b*x3\na == b", sep = "\n"),
                         lavaan::HolzingerSwineford1939, estimator = "MLM")
  uls <- lavaan::cfa(hs_model, lavaan::HolzingerSwineford1939,
                     estimator = "ULS", test = "scaled.shifted")
  for (fit in list(object, m0, m1_, m0_, general, uls)) {
    expect_true(uses_residual_factor(fit))
    df <- fit_df(fit)
    expect_equal(single_model_reduced_lambdas(fit, df), lavaan_lambdas(fit, df),
                 tolerance = 1e-8)
    if (is_classic_nt(fit)) {
      reference <- lapply(ugamma(fit, 3), ugamma_eigenvalues, df = df)
      expect_equal(single_model_reduced_lambdas(fit, df, 3), reference,
                   tolerance = 1e-8)
    }
  }
})

test_that("observed model information retains the lavaan reference route", {
  fit <- lavaan::cfa(hs_model, lavaan::HolzingerSwineford1939,
                     estimator = "MLR", information = "observed")
  expect_false(uses_residual_factor(fit))
  expect_equal(compute_spectra(fit, NULL, 1, fit_df(fit), "delta", "observed"),
               lavaan_lambdas(fit, fit_df(fit)))
})

test_that("FIML residual spectrum agrees with its full observed construction", {
  dat <- lavaan::HolzingerSwineford1939
  dat[seq(2, nrow(dat), 7), "x2"] <- NA
  for (group in list(NULL, "school")) {
    fit <- lavaan::cfa(hs_model, dat, group = group,
                       group.equal = if (is.null(group)) "" else "loadings",
                       estimator = "MLR", missing = "fiml")
    H <- fiml_h1_information_observed(fit)
    G <- fiml_gamma_matrix(fit)
    D <- model_effective_delta(fit)
    HD <- H %*% D
    U <- H - HD %*% generalized_inverse(crossprod(D, HD)) %*% t(HD)
    df <- fit_df(fit)
    expect_equal(fiml_lambdas(fit, df)$ug_biased,
                 fiml_sandwich_eigenvalues(U, G, df), tolerance = 1e-8)
  }
})
