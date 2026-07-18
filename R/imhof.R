# Tail probability of a positive (or mixed-sign) linear combination of
# independent chi-square_1 variables,
#
#     Q = sum_j lambda_j * Z_j^2,   Z_j ~ N(0, 1) independent,
#
# i.e. P(Q > q). This is the reference law for the eigenvalue-based SEM fit
# statistics: `lambdas` are the eigenvalues of the U*Gamma matrix and `q` is the
# observed chi-square. It is a drop-in replacement for
# `CompQuadForm::imhof(q, lambdas)$Qq`, so the package no longer depends on
# CompQuadForm.
#
# Three evaluations are combined, each chosen for the regime it is best in:
#
#   * Ruben (1962) series of central chi-squares -- used when every lambda_j is
#     positive (the usual case: `lambdas` are eigenvalues of U*Gamma). Exact to
#     machine precision and, crucially, it does NOT degrade in the tail, so the
#     small p-values produced by a poorly fitting model are accurate where
#     CompQuadForm is not. (A badly fitting textbook CFA already lands near
#     1e-7, well inside the unreliable zone.)
#   * Imhof (1961) numerical inversion of the characteristic function -- handles
#     mixed-sign eigenvalues (which do occur for some U*Gamma matrices, and are
#     why the statistics use `imhof` rather than the positive-only
#     `farebrother`). Exact up to quadrature error in the body of the
#     distribution; this is what CompQuadForm's `imhof` computes, and the two
#     agree to many digits there.
#   * Lugannani-Rice (1980) saddlepoint approximation -- the fallback once the
#     mixed-sign Imhof integral can no longer resolve the tail. CompQuadForm's
#     `imhof` returns quadrature noise (often *negative* probabilities) in that
#     regime and `davies` underflows to exactly 0; the saddlepoint instead keeps
#     a small *relative* error arbitrarily far into the tail.

# ---------------------------------------------------------------------------
# Imhof integral
# ---------------------------------------------------------------------------

# Integrand sin(theta(u)) / (u rho(u)) for h = 1, delta = 0, vectorised over u.
.imhof_integrand <- function(u, q, lambda) {
  lu <- outer(u, lambda)              # length(u) x r
  theta  <- 0.5 * rowSums(atan(lu)) - 0.5 * q * u
  logrho <- 0.25 * rowSums(log1p(lu^2))
  val <- sin(theta) * exp(-logrho) / u
  # limit as u -> 0: sin(theta)/u -> 0.5 (sum lambda - q)
  val[u == 0] <- 0.5 * (sum(lambda) - q)
  val
}

# Imhof upper tail at a single q. Returns list(prob, ok); `ok` is FALSE when the
# quadrature can no longer resolve the (tiny, oscillatory) tail.
.imhof_one <- function(q, lambda, rel.tol = 1e-10, subdivisions = 10000L) {
  res <- tryCatch(
    stats::integrate(.imhof_integrand, 0, Inf, q = q, lambda = lambda,
                     rel.tol = rel.tol, abs.tol = 1e-10,
                     subdivisions = subdivisions, stop.on.error = FALSE),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(prob = NA_real_, ok = FALSE))
  prob <- 0.5 + res$value / pi
  abserr <- res$abs.error / pi
  # Trust the QUADPACK error estimate rather than the convergence message: an
  # oscillatory integrand routinely returns an accurate value alongside a "max
  # subdivisions" message, while in the deep tail the reported error correctly
  # blows up relative to the (tiny) probability and we hand off to the
  # saddlepoint.
  ok <- is.finite(prob) && prob >= -1e-9 && prob <= 1 + 1e-9 &&
        abserr <= 1e-2 * max(prob, 1e-300)
  list(prob = min(max(prob, 0), 1), ok = ok)
}

# ---------------------------------------------------------------------------
# Lugannani-Rice saddlepoint (cumulant generating function of Q)
# ---------------------------------------------------------------------------

.K   <- function(t, lambda) sum(-0.5 * log(1 - 2 * lambda * t))
.Kp  <- function(t, lambda) sum(lambda / (1 - 2 * lambda * t))
.Kpp <- function(t, lambda) sum(2 * lambda^2 / (1 - 2 * lambda * t)^2)

# Solve K'(s) = q on the open domain where K is finite (K' is increasing).
.saddle_root <- function(q, lambda) {
  pos <- lambda[lambda > 0]
  neg <- lambda[lambda < 0]
  ub <- if (length(pos)) 1 / (2 * max(pos)) else Inf
  lb <- if (length(neg)) 1 / (2 * min(neg)) else -Inf
  Kp <- function(t) .Kp(t, lambda)
  span <- if (is.finite(ub) && is.finite(lb)) ub - lb else 1
  lo <- if (is.finite(lb)) lb + span * 1e-6 else -1
  hi <- if (is.finite(ub)) ub - span * 1e-6 else  1
  it <- 0L
  while (Kp(lo) > q && it < 300L) {
    lo <- if (is.finite(lb)) lb + (lo - lb) * 0.5 else lo * 2; it <- it + 1L
  }
  it <- 0L
  while (Kp(hi) < q && it < 300L) {
    hi <- if (is.finite(ub)) ub - (ub - hi) * 0.5 else hi * 2; it <- it + 1L
  }
  if (Kp(lo) > q || Kp(hi) < q) return(NA_real_)
  stats::uniroot(function(t) Kp(t) - q, lower = lo, upper = hi,
                 tol = .Machine$double.eps^0.6)$root
}

# Saddlepoint upper tail at a single q.
.saddle_one <- function(q, lambda) {
  mean_q <- sum(lambda)
  if (abs(q - mean_q) < 1e-9 * max(1, abs(mean_q))) {
    # Removable singularity at the mean: F(mean) = 1/2 + skewness / (6 sqrt(2 pi)).
    k2 <- sum(2 * lambda^2)
    k3 <- sum(8 * lambda^3)
    return(0.5 - k3 / (6 * sqrt(2 * pi) * k2^1.5))
  }
  s <- .saddle_root(q, lambda)
  if (is.na(s)) return(NA_real_)
  w <- sign(s) * sqrt(2 * (s * q - .K(s, lambda)))
  v <- s * sqrt(.Kpp(s, lambda))
  stats::pnorm(w, lower.tail = FALSE) + stats::dnorm(w) * (1 / v - 1 / w)
}

# ---------------------------------------------------------------------------
# Ruben (1962) series of central chi-squares (strictly positive weights)
# ---------------------------------------------------------------------------

# Exact to machine precision, and -- unlike the Imhof integral and Davies --
# it does NOT degrade in the tail: the survival is accumulated as a sum of
# non-negative terms, so there is no cancellation. This is the workhorse for the
# usual case where the U*Gamma eigenvalues are all positive (the regime where a
# poorly fitting model pushes the p-value far below what Imhof/Davies resolve).
# Returns list(prob, ok); ok = FALSE if the geometric series has not converged
# within `maxit` terms (a very wide eigenvalue spread), in which case the caller
# falls back to the Imhof/saddlepoint path.
.ruben_one <- function(q, lambda, eps = 1e-12, maxit = 5000L) {
  if (q <= 0) return(list(prob = 1, ok = TRUE))
  beta  <- min(lambda)
  m     <- length(lambda)                 # total d.f.; each component is chi^2_1
  gamma <- 1 - beta / lambda              # in [0, 1)
  g     <- function(i) 0.5 * sum(gamma^(i + 1))
  acoef <- exp(0.5 * sum(log(beta / lambda)))             # a_0
  gcoef <- g(0L)
  surv  <- acoef * stats::pchisq(q / beta, df = m, lower.tail = FALSE)
  asum  <- acoef
  k     <- 1L
  repeat {
    # a_k = (1/k) sum_{i=0}^{k-1} g_i a_{k-1-i}
    ak    <- sum(gcoef * rev(acoef)) / k
    acoef <- c(acoef, ak)
    gcoef <- c(gcoef, g(k))
    surv  <- surv + ak * stats::pchisq(q / beta, df = m + 2L * k,
                                       lower.tail = FALSE)
    asum  <- asum + ak
    # The omitted weight 1 - asum bounds the remaining (non-negative) tail of
    # the series, so this is a genuine *relative* error guarantee on `surv`.
    if (1 - asum < eps * max(surv, .Machine$double.xmin))
      return(list(prob = min(max(surv, 0), 1), ok = TRUE))
    if (k >= maxit)
      return(list(prob = min(max(surv, 0), 1), ok = FALSE))
    k <- k + 1L
  }
}

# ---------------------------------------------------------------------------
# Front-end
# ---------------------------------------------------------------------------

#' Upper tail of a linear combination of chi-square_1 variables
#'
#' Computes \eqn{P(Q > q)} for \eqn{Q = \sum_j \lambda_j Z_j^2} with
#' independent \eqn{Z_j \sim N(0,1)}. Drop-in replacement for
#' `CompQuadForm::imhof(q, lambda)$Qq`: the Imhof integral in the body of the
#' distribution, a Lugannani-Rice saddlepoint in the far tail where the integral
#' degrades.
#'
#' @param q Numeric scalar (or vector) of thresholds; the observed chi-square.
#' @param lambda Numeric vector of eigenvalues (may be mixed sign).
#' @return Numeric vector of upper-tail probabilities, length `length(q)`.
#' @references
#'
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in normal
#' variables. *Biometrika*, 48(3/4), 419--426.
#' \doi{10.1093/biomet/48.3-4.419}
#'
#' Lugannani, R., & Rice, S. O. (1980). Saddle point approximation for the
#' distribution of the sum of independent random variables. *Advances in
#' Applied Probability*, 12(2), 475--490. \doi{10.2307/1426607}
#'
#' Ruben, H. (1962). Probability content of regions under spherical normal
#' distributions, IV: The distribution of homogeneous and non-homogeneous
#' quadratic functions of normal variables. *The Annals of Mathematical
#' Statistics*, 33(2), 542--570. \doi{10.1214/aoms/1177704580}
#' @keywords internal
imhof_pvalue <- function(q, lambda) {
  lambda <- lambda[is.finite(lambda)]
  if (length(lambda) == 0L)
    return(as.numeric(q < 0))                       # Q is degenerate at 0
  # Drop numerically negligible eigenvalues. They contribute nothing to Q, and
  # this keeps a noise-level negative eigenvalue (~ -1e-12 out of an O(1)
  # spectrum) from bumping the whole computation off the exact, positive-weight
  # Ruben path onto the mixed-sign one.
  lambda <- lambda[abs(lambda) > 1e-8 * max(abs(lambda))]
  if (length(lambda) == 0L)
    return(as.numeric(q < 0))
  if (length(lambda) == 1L)                         # Q = lambda * chi^2_1, exact
    return(stats::pchisq(q / lambda, df = 1,
                         lower.tail = lambda < 0))

  positive <- all(lambda > 0)
  one <- function(qi) {
    # Positive weights: Ruben's series is exact in every regime, including the
    # deep tail where a badly fitting model lands.
    if (positive) {
      rb <- .ruben_one(qi, lambda)
      if (rb$ok) return(rb$prob)
    }
    # Mixed-sign weights (or a non-converged Ruben series): the exact Imhof
    # integral in the body, the saddlepoint once the integral can no longer
    # resolve the tail.
    im <- .imhof_one(qi, lambda)
    if (im$ok) return(im$prob)
    sp <- .saddle_one(qi, lambda)
    if (is.finite(sp)) {
      return(min(max(sp, 0), 1))
    }
    im$prob
  }
  vapply(q, one, numeric(1))
}
