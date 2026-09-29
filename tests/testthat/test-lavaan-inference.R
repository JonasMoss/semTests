test_that("Gaussian score rows are derivatives of the actual log likelihood", {
  Sigma <- matrix(c(2, .4, .4, 1.5), 2)
  mu <- c(.7, -.4)
  X <- rbind(c(-1, .2), c(2, 1))
  pack <- gaussian_moment_geometry(Sigma, TRUE)
  theta <- c(mu, Sigma[pack$index])
  loglik <- function(par, x) {
    S <- matrix(0, 2, 2)
    S[pack$index] <- par[-(1:2)]
    S[upper.tri(S)] <- t(S)[upper.tri(S)]
    a <- x - par[1:2]
    -.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) -
      .5 * drop(crossprod(a, solve(S, a)))
  }
  actual <- gaussian_score_rows(sweep(X, 2, mu, "-"), pack)
  h <- 1e-6
  expected <- t(vapply(seq_len(nrow(X)), function(i) {
    vapply(seq_along(theta), function(j) {
      step <- numeric(length(theta)); step[j] <- h
      (loglik(theta + step, X[i, ]) - loglik(theta - step, X[i, ])) / (2 * h)
    }, numeric(1))
  }, numeric(length(theta))))
  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("Gaussian Fisher information matches a population likelihood Hessian", {
  S0 <- matrix(c(2, .4, .4, 1.5), 2)
  mu0 <- c(.7, -.4)
  pack <- gaussian_moment_geometry(S0, TRUE)
  theta <- c(mu0, S0[pack$index])
  criterion <- function(par) {
    S <- matrix(0, 2, 2); S[pack$index] <- par[-(1:2)]
    S[upper.tri(S)] <- t(S)[upper.tri(S)]
    a <- mu0 - par[1:2]
    .5 * (as.numeric(determinant(S, logarithm = TRUE)$modulus) +
      sum(diag(solve(S, S0 + tcrossprod(a)))))
  }
  h <- 1e-4
  H <- matrix(0, length(theta), length(theta))
  for (i in seq_along(theta)) for (j in seq_along(theta)) {
    ei <- ej <- numeric(length(theta)); ei[i] <- h; ej[j] <- h
    H[i, j] <- (criterion(theta + ei + ej) - criterion(theta + ei - ej) -
      criterion(theta - ei + ej) + criterion(theta - ei - ej)) / (4 * h^2)
  }
  expect_equal(pack$information, H, tolerance = 1e-6)
})

test_that("global score projection removes the model tangent across groups", {
  for (fit in list(object, m0, m0_no_groups)) {
    out <- score_test(fit, keep = TRUE)
    G <- out$ingredients$projection
    H <- lavaan::lav_matrix_bdiag(out$ingredients$information_blocks)
    D <- out$ingredients$tangent
    expect_equal(unname(crossprod(G, H %*% G)), diag(out$df), tolerance = 1e-9)
    expect_lt(max(abs(crossprod(D, H %*% G))), 1e-9)
    expect_equal(out$df, as.integer(fit_df(fit)))
    expect_lt(out$diagnostics$nuisance_stationarity_norm, 1e-4)
    parts <- ml_global_score_parts(fit, "uncentered", tile_rows = 1L)
    expect_equal(parts$score, out$ingredients$score, tolerance = 1e-9)
    expect_equal(parts$score_covariance, out$ingredients$score_covariance, tolerance = 1e-9)
    centered <- score_test(fit, score_covariance = "centered", keep = TRUE)
    expect_equal(centered$statistic, out$statistic)
    expect_gte(min(eigen(out$ingredients$score_covariance - centered$ingredients$score_covariance,
                         symmetric = TRUE, only.values = TRUE)$values), -1e-10)
    expect_identical(score_test(fit, metric = "sandwich")$reference, "chisq")
  }
})

test_that("adapters reject unsupported inference instead of changing its meaning", {
  expect_error(score_test(m1_), "complete-data normal-likelihood ML")
  fit <- object
  fit@Data@std.ov <- TRUE
  expect_error(score_test(fit), "complete-data normal-likelihood ML")
  fit <- object
  fit@Model@ceq.nonlinear.idx <- 1L
  expect_error(score_test(fit), "affine equality")
})
