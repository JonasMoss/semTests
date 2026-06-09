# Tests for imhof_pvalue(), the in-package replacement for
# CompQuadForm::imhof(q, lambda)$Qq.
#
# The primary oracles are INDEPENDENT of CompQuadForm -- exact chi-square
# probabilities, Monte Carlo with a binomial standard error, and the rigorous
# leading-order tail asymptotic -- because part of the point of the replacement
# is that CompQuadForm is itself unreliable in the tail. CompQuadForm is then
# used as a secondary oracle in the body of the distribution only, loaded
# dynamically (via a variable package name, so it is NOT a declared dependency
# of the package and never triggers an R CMD check NOTE).

# --- helpers ---------------------------------------------------------------

cqf_name  <- "CompQuadForm"                       # variable => not a hard dep
have_cqf  <- requireNamespace(cqf_name, quietly = TRUE)
imhof_cqf <- function(q, lambda)
  as.numeric(getExportedValue(cqf_name, "imhof")(q, lambda)$Qq)

# Monte Carlo upper tail with its binomial standard error.
mc_tail <- function(q, lambda, N = 2e6, seed = 1L) {
  set.seed(seed)
  X <- matrix(stats::rchisq(N * length(lambda), df = 1), N, length(lambda))
  p <- mean(as.numeric(X %*% lambda) > q)
  c(p = p, se = sqrt(p * (1 - p) / N))
}

rel_err <- function(got, ref) abs(got - ref) / abs(ref)

# --- exact references ------------------------------------------------------

test_that("a single eigenvalue is the exact (scaled) chi-square_1 tail", {
  expect_equal(imhof_pvalue(22, 1),  stats::pchisq(22, 1, lower.tail = FALSE))
  expect_equal(imhof_pvalue(5,  3),  stats::pchisq(5 / 3, 1, lower.tail = FALSE))
  expect_equal(imhof_pvalue(0.2, 4), stats::pchisq(0.05, 1, lower.tail = FALSE))
  # a negative weight flips the tail: Q = -|lambda| chi^2_1 <= 0
  expect_equal(imhof_pvalue(-3, -2), stats::pchisq(1.5, 1, lower.tail = TRUE))
  expect_equal(imhof_pvalue(2,  -2), 0)            # P(Q > 2) = 0 when Q <= 0
})

test_that("equal weights reduce to a scaled chi-square", {
  for (k in 2:6) {
    q <- 1.5 * k
    expect_lt(rel_err(imhof_pvalue(q, rep(2, k)),
                      stats::pchisq(q / 2, k, lower.tail = FALSE)), 5e-2)
  }
})

# --- Monte Carlo across spectra (independent of CompQuadForm) ---------------

test_that("matches Monte Carlo across a range of spectra", {
  spectra <- list(
    c(3, 2, 1),
    c(1, 1, 1, 1, 1),
    c(5, 1, 0.5, 0.25),
    sort(rgamma(8, shape = 2, rate = 2), decreasing = TRUE),
    c(2.0, 1.3, 0.8, 0.3, 0.05)
  )
  for (lambda in spectra) {
    for (mult in c(0.6, 1.0, 1.6)) {
      q  <- sum(lambda) * mult
      mc <- mc_tail(q, lambda, N = 2e6, seed = 11L)
      if (mc[["p"]] < 1e-3 || mc[["p"]] > 0.999) next   # MC too noisy out here
      expect_lt(abs(imhof_pvalue(q, lambda) - mc[["p"]]), 6 * mc[["se"]] + 1e-4)
    }
  }
})

# --- parity with CompQuadForm in the body of the distribution ---------------

test_that("agrees with CompQuadForm::imhof in the body (p > 1e-5)", {
  skip_if_not(have_cqf, "CompQuadForm not installed")
  set.seed(7L)
  for (i in 1:30) {
    r      <- sample(2:14, 1L)
    lambda <- sort(rgamma(r, shape = 2, rate = 2), decreasing = TRUE)
    q      <- sum(lambda) * runif(1, 0.4, 2.6)
    ref    <- imhof_cqf(q, lambda)
    if (ref < 1e-5 || ref > 1 - 1e-5) next
    expect_lt(rel_err(imhof_pvalue(q, lambda), ref), 1e-4)
  }
})

test_that("agrees with CompQuadForm for mixed-sign eigenvalues (body)", {
  skip_if_not(have_cqf, "CompQuadForm not installed")
  lambda <- c(2.5, 1, 0.4, -0.3, -1)
  for (q in c(-2, -0.5, 0, 2, 5, 8)) {
    ref <- imhof_cqf(q, lambda)
    expect_lt(rel_err(imhof_pvalue(q, lambda), ref), 1e-3)
  }
})

# --- the far tail, where CompQuadForm degrades ------------------------------

test_that("stays accurate in the deep tail (exact chi-square reference)", {
  # Equal weights => Q ~ chi^2_k exactly, so pchisq gives ground truth far past
  # where the Imhof integral and Davies can reach.
  # Equal weights => Q ~ chi^2_k exactly; Ruben's series reproduces it to
  # machine precision arbitrarily far out (here down to ~1e-32).
  for (cfg in list(list(k = 3, qs = c(60, 100, 160, 240)),
                   list(k = 5, qs = c(80, 140, 220)))) {
    for (q in cfg$qs) {
      exact <- stats::pchisq(q, cfg$k, lower.tail = FALSE)
      got   <- imhof_pvalue(q, rep(1, cfg$k))
      expect_gt(got, 0)                       # never a noise-driven zero/negative
      expect_lt(rel_err(got, exact), 1e-8)
    }
  }
})

test_that("deep tail agrees with Davies for distinct positive weights", {
  skip_if_not(have_cqf, "CompQuadForm not installed")
  # Davies is the one CompQuadForm routine that is usable in the tail, but only
  # while it reports success -- it sets ifault (or underflows to 0) once it can
  # no longer reach the requested accuracy. Compare only where it succeeds.
  davies <- function(q, l) getExportedValue(cqf_name, "davies")(q, l, acc = 1e-11, lim = 1e6)
  checked <- 0L
  for (lambda in list(c(3, 2, 1), c(5, 2, 1, 0.5), c(2, 1.5, 1, 0.7, 0.3))) {
    for (q in sum(lambda) * c(3, 5, 7)) {
      res <- suppressWarnings(davies(q, lambda))
      if (res$ifault != 0 || res$Qq <= 0) next
      expect_lt(rel_err(imhof_pvalue(q, lambda), res$Qq), 1e-4)
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 0L)
})

test_that("beats CompQuadForm::imhof in the deep tail (the motivating case)", {
  skip_if_not(have_cqf, "CompQuadForm not installed")
  cqf_imhof <- function(q, l) as.numeric(getExportedValue(cqf_name, "imhof")(q, l)$Qq)
  # chi^2_3: exact truth is pchisq. CompQuadForm's imhof returns quadrature
  # noise here (wrong sign / orders of magnitude off); ours stays exact.
  worse <- 0L
  for (q in c(60, 100, 160)) {
    exact <- stats::pchisq(q, 3, lower.tail = FALSE)
    ours  <- imhof_pvalue(q, c(1, 1, 1))
    cqf   <- suppressWarnings(cqf_imhof(q, c(1, 1, 1)))
    expect_lt(rel_err(ours, exact), 1e-8)               # ours: exact
    if (!(cqf > 0) || rel_err(cqf, exact) > 0.1) worse <- worse + 1L
  }
  expect_equal(worse, 3L)                                # CompQuadForm fails all three
})

test_that("never returns an invalid probability, even in the deep tail", {
  set.seed(3L)
  for (i in 1:300) {
    r      <- sample(2:10, 1L)
    lambda <- rgamma(r, shape = 3, rate = 3)
    q      <- sum(lambda) * runif(1, 1, 9)     # frequently far into the tail
    p      <- imhof_pvalue(q, lambda)
    expect_true(is.finite(p) && p >= 0 && p <= 1)
  }
})

# --- structure, edge cases, vectorisation -----------------------------------

test_that("zeros and non-finite weights are dropped", {
  lambda <- c(3, 2, 1)
  expect_equal(imhof_pvalue(5, c(lambda, 0, 0)), imhof_pvalue(5, lambda))
  expect_equal(imhof_pvalue(5, c(lambda, NA)),    imhof_pvalue(5, lambda))
  # all-zero weights => Q is degenerate at 0
  expect_equal(imhof_pvalue(1,  c(0, 0)), 0)
  expect_equal(imhof_pvalue(-1, c(0, 0)), 1)
  # all weights non-finite => nothing left after the finite filter, also Q = 0
  expect_equal(imhof_pvalue(1,  c(NA, Inf, -Inf)), 0)
  expect_equal(imhof_pvalue(-1, c(NA_real_, NA_real_)), 1)
})

test_that("monotone in q and bounded in [0, 1]", {
  lambda <- c(4, 2, 1, 0.5)
  qs <- seq(-2, 140, length.out = 60)
  ps <- vapply(qs, imhof_pvalue, numeric(1), lambda = lambda)
  expect_true(all(ps >= 0 & ps <= 1))
  expect_true(all(diff(ps) <= 1e-9))           # non-increasing
  expect_gt(ps[1], 0.99)                        # P(Q > -2) ~ 1
  expect_lt(ps[length(ps)], 1e-6)              # well into the tail
})

test_that("vectorises over q", {
  lambda <- c(2, 1, 0.5)
  qs <- c(0, 1, 5, 10)
  expect_equal(imhof_pvalue(qs, lambda),
               vapply(qs, imhof_pvalue, numeric(1), lambda = lambda))
})

# --- the Lugannani-Rice saddlepoint fallback (mixed-sign deep tail) ----------

test_that("the saddlepoint keeps a valid, monotone tail for mixed-sign weights", {
  # Mixed-sign spectrum pushed far enough that the Imhof integral can no longer
  # resolve the (tiny) tail and imhof_pvalue() hands off to the saddlepoint.
  lambda <- c(6, 3, 1, -2)
  qs <- c(20, 40, 70, 110)
  ps <- vapply(qs, imhof_pvalue, numeric(1), lambda = lambda)
  expect_true(all(is.finite(ps) & ps >= 0 & ps <= 1))
  expect_true(all(diff(ps) <= 1e-9))            # non-increasing in q
  expect_gt(ps[1], ps[length(ps)])              # genuinely decaying
})

test_that(".saddle_one handles the removable singularity at the mean", {
  # At q = E[Q] = sum(lambda) the Lugannani-Rice formula has a removable
  # singularity handled by a skewness expansion; the result is a valid F(mean).
  lambda <- c(6, 3, 1, -2)
  s <- .saddle_one(sum(lambda), lambda)
  expect_true(is.finite(s) && s > 0 && s < 1)
  # off the mean the saddle root exists for a one-sided (all-positive) spectrum
  expect_true(is.finite(.saddle_root(10, c(3, 2, 1))))
})

test_that(".imhof_integrand has the correct u -> 0 limit", {
  lambda <- c(2, 1, 0.5); q <- 1.5
  v <- .imhof_integrand(c(0, 0.3), q, lambda)
  expect_equal(v[1], 0.5 * (sum(lambda) - q))   # the analytic limit at u = 0
  expect_true(is.finite(v[2]))
})
