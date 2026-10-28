# Refusals, guards, and fallbacks that ordinary use rarely reaches. Most guard
# against inconsistent inputs or a changed lavaan; a few are reached only by
# altering a fitted object, which stands in for a fit lavaan would not produce.

hs <- lavaan::HolzingerSwineford1939
one_factor <- "f =~ x1 + a*x2 + b*x3 + x4"
three_factor <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"

test_that("inference_control() validates and prints its settings", {
  invalid <- "semTests_error_invalid_inference"
  expect_error(inference_control(derivative_step = 1), "derivative_step", class = invalid)
  expect_error(inference_control(tolerance = 2), "tolerance", class = invalid)
  expect_error(inference_control(keep = NA), "keep", class = invalid)
  expect_output(print(inference_control()), "curvature +<by test>")
  expect_output(print(inference_control(curvature = "observed")), "curvature +observed")

  fit <- lavaan::cfa(one_factor, hs)
  expect_error(semtest(fit, control = list()), "inference_control", class = invalid)
  expect_error(semtest(fit, unknown = 1), "setting: unknown\\.", class = invalid)
  expect_error(semtest(fit, unknown = 1, other = 2), "settings: unknown, other",
               class = invalid)
})

test_that("semtest() explains when the criterion statistic is unavailable", {
  invalid <- "semTests_error_invalid_inference"
  fit <- lavaan::cfa(one_factor, hs)
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), hs)
  template <- lavaan::cfa(one_factor, hs, do.fit = FALSE)
  expect_error(semtest(fit, statistic = "criterion"), "use pvalues", class = invalid)
  expect_error(semtest(restricted, template, statistic = "criterion"),
               "must be optimized", class = invalid)
})

test_that("overall score tests refuse fits without a testable complement", {
  saturated <- lavaan::cfa("f =~ x1 + x2 + x3", hs)
  expect_error(semtest(saturated), "no testable saturated complement",
               class = "semTests_error_invalid_inference")
  grouped <- lavaan::cfa("f =~ x1 + x2 + x3 + x4", hs, group = "school",
                         group.w.free = TRUE)
  expect_error(semtest(grouped), "requires complete-data",
               class = "semTests_error_unsupported_inference")
})

test_that("overall score tests check the data against the fitted moments", {
  fit <- lavaan::cfa(one_factor, hs)
  altered <- fit
  altered@SampleStats@nobs[[1]] <- 300L
  expect_error(semtest(altered), "ordering do not agree",
               class = "semTests_error_invalid_inference")
  expect_error(gaussian_pattern_parts(altered), "case counts do not agree")
  local_mocked_bindings(residual_metric_factor = function(...) NULL)
  expect_error(semtest(fit), "positive definite",
               class = "semTests_error_unstable_spectrum")
})

test_that("criterion inference refuses unsupported least-squares fits", {
  sample_cov <- stats::cov(hs[, c("x1", "x2", "x3", "x4")])
  fit <- lavaan::cfa(one_factor, sample.cov = sample_cov, sample.nobs = 301,
                     estimator = "GLS")
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"),
                            sample.cov = sample_cov, sample.nobs = 301,
                            estimator = "GLS")
  expect_error(semtest(restricted, fit), "raw-data ML or complete-data LS")

  bounded <- lavaan::cfa(paste(one_factor, "a > 0", sep = "\n"), hs, estimator = "GLS")
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), hs, estimator = "GLS")
  expect_error(semtest(bounded, restricted), "affine equality constraints only")

  mixed <- hs[, c("x1", "x2", "x3", "x4")]
  mixed[c("x1", "x2")] <- lapply(mixed[c("x1", "x2")], function(x) as.integer(x > median(x)))
  fit <- lavaan::cfa(one_factor, mixed, ordered = c("x1", "x2"))
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), mixed,
                            ordered = c("x1", "x2"))
  expect_error(semtest(restricted, fit), "requires all variables ordered")
})

test_that("nested score tests check the embedding of the restricted model", {
  fit <- lavaan::cfa(one_factor, hs)
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), hs)

  gls <- lavaan::cfa(one_factor, hs, estimator = "GLS", do.fit = FALSE)
  expect_error(semtest(restricted, gls), "same fitting criterion")

  bounded <- lavaan::cfa(paste(one_factor, "a > 0", sep = "\n"), hs)
  expect_error(semtest(restricted, bounded), "only affine equality constraints")

  wider <- lavaan::cfa(paste(one_factor, "x1 ~~ x2", sep = "\n"), hs, do.fit = FALSE)
  expect_error(semtest(restricted, wider), "every alternative parameter")

  ambiguous <- fit
  ambiguous@ParTable$free[ambiguous@ParTable$free == 2L] <- 1L
  expect_error(semtest(restricted, ambiguous), "ordering is ambiguous")

  offset <- lavaan::cfa(paste(one_factor, "a == b + .1", sep = "\n"), hs,
                        do.fit = FALSE)
  expect_error(semtest(restricted, offset), "does not satisfy alternative affine")

  calls <- 0L
  real_implied <- inference_implied
  local_mocked_bindings(inference_implied = function(object) {
    calls <<- calls + 1L
    out <- real_implied(object)
    if (calls == 2L) out[[1]] <- out[[1]] + 1
    out
  })
  expect_error(semtest(restricted, fit), "parameterizations do not align")
})

test_that("nested score tests check the moment tangents", {
  fit <- lavaan::cfa(one_factor, hs)
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), hs)
  calls <- 0L
  real_delta <- inference_delta
  local_mocked_bindings(inference_delta = function(object) {
    calls <<- calls + 1L
    out <- real_delta(object)
    if (calls == 2L) out[[1]] <- out[[1]] * 2
    out
  })
  expect_error(semtest(restricted, fit), "tangents do not align")
})

test_that("categorical nested score tests need a shared parameterization", {
  d <- hs[, c("x1", "x2", "x3", "x4")]
  d[] <- lapply(d, function(x) as.integer(x > median(x)))
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), d,
                            ordered = names(d))
  theta <- lavaan::cfa(one_factor, d, ordered = names(d),
                       parameterization = "theta", do.fit = FALSE)
  expect_error(semtest(restricted, theta), "categorical parameterization")
})

test_that("criterion differences check the first stage and the optimization", {
  fit <- lavaan::cfa(one_factor, hs, estimator = "GLS")
  restricted <- lavaan::cfa(paste(one_factor, "a == b", sep = "\n"), hs, estimator = "GLS")
  altered <- fit
  altered@SampleStats@WLS.V[[1]] <- altered@SampleStats@WLS.V[[1]] * 1.01
  expect_error(semtest(restricted, altered, statistic = "criterion"),
               "identical first-stage moments")
  unoptimized <- restricted
  unoptimized@optim$fx <- fit@optim$fx - 1
  expect_error(semtest(unoptimized, fit, statistic = "criterion"),
               "smaller than the alternative")
})

test_that("score matrices refuse unusable inputs", {
  expect_identical(inference_complement(matrix(0, 3, 0)), diag(3))
  expect_error(inference_complement(cbind(c(1, 0, 0), 0)), "zero columns")
  expect_error(profile_score_test_matrices(c(1, 1), diag(2), -diag(2), matrix(0, 2, 0)),
               "positive semidefinite")
  expect_error(profile_score_test_matrices(c(1, 1), diag(2), diag(2), c(1, 0)),
               "finite conformable matrix")
  expect_error(inference_chol(matrix(NA_real_, 1, 1), "x"), "non-finite",
               class = "semTests_error_invalid_inference")
  expect_error(score_test_matrices(c(1e200, 1e200), diag(c(1e-200, 1e-200)), diag(2)),
               "overflowed", class = "semTests_error_invalid_inference")
  expect_error(score_test_matrices(c(1, 1), diag(c(1, 1e300)), diag(c(1, 1e-30))),
               "numerically singular", class = "semTests_error_unstable_spectrum")
})

test_that("numerical criterion curvature checks its step and symmetry", {
  expect_error(criterion_hessian(identity, 1, matrix(1), 1), "derivative_step")
  skewed <- function(theta) c(theta[1] + 2 * theta[2], theta[2])
  expect_error(criterion_hessian(skewed, c(1, 1), diag(2), 1e-4), "asymmetric")
})

test_that("residual factorizations check their dimensions", {
  expect_error(residual_metric_factor(list(diag(2)), matrix(1, 3, 1), 1),
               "dimensions do not agree")
  expect_error(residual_metric_factor(list(matrix(1, 2, 3)), matrix(1, 2, 1), 1),
               "finite square matrix")
  expect_error(residual_metric_factor(list(diag(2)), matrix(c(NA, 1), 2, 1), 1),
               "non-finite")
  expect_error(residual_metric_factor(list(diag(3)), cbind(1:3, 0), 1),
               "rank deficient")
  expect_error(positive_symmetric_eigenvalues(diag(2), 3, "test"), "invalid dimensions")
  # A singular C keeps the general construction.
  expect_equal(restriction_spectrum(diag(c(1, 0)), diag(2), 1, "test"), 1)
})

test_that("full-matrix spectra agree with the residual factorization", {
  fit <- lavaan::cfa(three_factor, hs)
  missing <- hs
  set.seed(1)
  missing$x1[sample(nrow(missing), 40)] <- NA
  fiml <- lavaan::cfa(three_factor, missing, missing = "fiml")
  factored <- list(pvalues(fit, "PEBA4_ML"), pvalues(fit, "SB_UG_ML"), pvalues(fiml))
  local_mocked_bindings(residual_metric_factor = function(...) NULL)
  expect_equal(pvalues(fit, "PEBA4_ML"), factored[[1]], tolerance = 1e-8)
  expect_equal(pvalues(fit, "SB_UG_ML"), factored[[2]], tolerance = 1e-8)
  expect_equal(pvalues(fiml), factored[[3]], tolerance = 1e-8)
})

test_that("observed-information ML fits use the unbiased gamma directly", {
  fit <- lavaan::cfa(three_factor, hs, information = "observed")
  expected <- lavaan::cfa(three_factor, hs)
  expect_equal(as.numeric(pvalues(fit, "SB_UG_ML")),
               as.numeric(pvalues(expected, "SB_UG_ML")), tolerance = 1e-2)
})

test_that("internal gamma and spectrum helpers refuse the unbiased gamma off ML", {
  gls <- lavaan::cfa(three_factor, hs, estimator = "GLS")
  expect_error(gamma_matrices(gls, 2), "only defined for continuous")
  expect_null(gamma_rescale(list(), gls))
  missing <- hs
  set.seed(2)
  missing$x1[sample(nrow(missing), 40)] <- NA
  fiml <- lavaan::cfa(three_factor, missing, missing = "fiml")
  expect_error(compute_spectra(fiml, NULL, 2, 24L, "delta", "observed"),
               "not defined for FIML")
  expect_error(compute_spectra(fiml, fiml, 2, 1L, "delta", "observed"),
               "not defined for FIML")
})

test_that("printing keeps nested results saved without A.method", {
  m1 <- lavaan::cfa(three_factor, hs)
  m0 <- lavaan::cfa(sub("x2 + x3", "a*x2 + a*x3", three_factor, fixed = TRUE), hs)
  result <- pvalues_nested(m0, m1)
  info <- attr(result, "semtests")
  info$A.method <- NULL
  attr(result, "semtests") <- info
  expect_output(print(result), "nested \\(method 2000\\)$")
})

test_that("moment projection checks its blocks against the fit", {
  fit <- lavaan::cfa(three_factor, hs)
  expect_error(projected_gamma(fit, matrix(1, 3, 1)), "dimensions do not agree")
  altered <- fit
  altered@SampleStats@nobs[[1]] <- 300L
  F <- matrix(1, 45, 1)
  expect_error(projected_gamma(altered, F), "sample sizes do not agree")
})

test_that("ordinal pair scores handle infinite edges and refuse empty cells", {
  expect_identical(bvn_density(c(Inf, -Inf), c(1, Inf), .3), c(0, 0))
  out <- pair_scores(c(1L, 1L), c(1L, 2L), 0, 0, .2)
  expect_true(all(is.finite(unlist(out))))
  expect_error(pair_scores(2L, 1L, c(0, 0), 0, .2), "nonempty model-implied")
})

test_that("score tests refuse non-finite raw observations", {
  fit <- lavaan::cfa(one_factor, hs)
  altered <- fit
  altered@Data@X[[1]][1, 1] <- Inf
  expect_error(gaussian_pattern_parts(altered), "ordering do not agree")
})

test_that("least-squares first stages check their inputs", {
  fit <- lavaan::cfa(one_factor, hs, estimator = "GLS")
  altered <- fit
  altered@Data@X[[1]] <- altered@Data@X[[1]] * 1.01
  expect_error(continuous_first_stage(altered, 1, "estimated"),
               "does not match the supported empirical weight map")
  expect_error(ls_first_stage(fit, "estimated", 1), "derivative_step")
  rescaled <- lavaan::cfa(one_factor, hs, estimator = "GLS", gamma.n.minus.one = TRUE)
  expect_error(ls_first_stage(rescaled, "estimated", 1e-4),
               "alternative Gamma normalization")

  d <- hs[, c("x1", "x2", "x3", "x4")]
  d[] <- lapply(d, function(x) as.integer(cut(x, quantile(x, c(0, 1/3, 2/3, 1)),
                                               include.lowest = TRUE)))
  uls <- lavaan::cfa(one_factor, d, ordered = names(d), estimator = "ULS")
  first <- ordinal_first_stage(uls, 1, "estimated", 1e-4)
  expect_equal(first$correction(first$observed), matrix(0, nrow(d), length(first$observed)))
  incomplete <- uls
  incomplete@Data@X[[1]][1, 1] <- NA
  expect_error(ordinal_first_stage(incomplete, 1, "fixed", 1e-4), "complete observations")
  gap <- uls
  X <- gap@Data@X[[1]]
  X[X[, 1] == 2, 1] <- 1
  X[X[, 1] == 3, 1] <- 3
  gap@Data@X[[1]] <- X
  expect_error(ordinal_first_stage(gap, 1, "fixed", 1e-4), "consecutive and nonempty")
})
