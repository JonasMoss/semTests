# Analytic influence function of the ordinal first stage.
#
# The first stage estimates thresholds and polychoric correlations in two
# conditional steps: thresholds from each variable's marginal multinomial, then
# each polychoric correlation by ML at those fixed thresholds. Stacking the two
# gives a just-identified estimating equation whose influence function is
#
#   phi_i = n * SC_i %*% t(B_inv),   B_inv = A^{-1}
#
# with A block lower triangular because the thresholds do not depend on the
# correlations. This reproduces the convention lavaan uses for its categorical
# Gamma, so the resulting weight matrices agree with lavaan's.
#
# The alternative is to differentiate the re-solved first-stage map numerically.
# That targets the exact Jacobian rather than the outer-product bread, agrees to
# O(n^-1/2), and costs a full polychoric refit per response pattern.

# Bivariate normal density. Zero at infinite arguments, which is what the
# rectangle corners need.
bvn_density <- function(h, k, rho) {
  out <- numeric(length(h))
  ok <- is.finite(h) & is.finite(k)
  if (!any(ok)) return(out)
  hh <- h[ok]
  kk <- k[ok]
  d <- 1 - rho^2
  out[ok] <- exp(-(hh^2 - 2 * rho * hh * kk + kk^2) / (2 * d)) /
    (2 * pi * sqrt(d))
  out
}

# Bivariate normal CDF with explicit handling of infinite arguments, which
# pbivnorm does not accept.
bvn_cdf <- function(h, k, rho) {
  out <- numeric(length(h))
  neg <- (h == -Inf) | (k == -Inf)
  out[neg] <- 0
  both <- (h == Inf) & (k == Inf)
  out[both] <- 1
  hi <- (h == Inf) & !neg & !both
  out[hi] <- stats::pnorm(k[hi])
  ki <- (k == Inf) & !neg & !both & !hi
  out[ki] <- stats::pnorm(h[ki])
  rest <- !neg & !both & !hi & !ki
  if (any(rest)) {
    out[rest] <- pbivnorm::pbivnorm(h[rest], k[rest], rho = rho)
  }
  out
}

# Probability of a threshold rectangle, and its correlation derivative. The
# derivative uses Plackett's identity dPhi2/drho = phi2, so no quadrature.
bvn_rect <- function(h_lo, h_hi, k_lo, k_hi, rho, derivative = FALSE) {
  f <- if (derivative) bvn_density else bvn_cdf
  f(h_hi, k_hi, rho) - f(h_lo, k_hi, rho) -
    f(h_hi, k_lo, rho) + f(h_lo, k_lo, rho)
}

# Per-case influence of one variable's thresholds. tau_a = qnorm(cumulative
# proportion), so the delta method gives this directly and exactly; there is no
# bread to invert. Using the multinomial outer product here instead would be a
# different and wrong object, because the categories share thresholds.
threshold_influence <- function(codes, tau) {
  vapply(seq_along(tau), function(a) {
    ((codes <= a) - stats::pnorm(tau[a])) / stats::dnorm(tau[a])
  }, numeric(length(codes)))
}

# Per-case scores of one polychoric correlation, and of that pair's cell
# probability with respect to each variable's thresholds.
#
# The threshold channel is NOT the marginal multinomial score: it is the pair
# likelihood's own derivative. Only this one satisfies the information identity
# against the correlation score, which is what makes the outer product a valid
# estimate of the cross bread block.
#
#   dP/dtau = phi(tau) * [Phi(z(hi)) - Phi(z(lo))],  z(t) = (t - rho*tau)/sqrt(1-rho^2)
#
# evaluated at the cell's own edges, so each case touches at most two thresholds
# per variable: the upper edge of its category and the lower one.
pair_scores <- function(code_u, code_v, tau_u, tau_v, rho) {
  bounds_u <- c(-Inf, tau_u, Inf)
  bounds_v <- c(-Inf, tau_v, Inf)
  lo_u <- bounds_u[code_u]
  hi_u <- bounds_u[code_u + 1L]
  lo_v <- bounds_v[code_v]
  hi_v <- bounds_v[code_v + 1L]
  probability <- bvn_rect(lo_u, hi_u, lo_v, hi_v, rho)
  if (any(probability <= 0)) {
    stop("Ordinal influence requires nonempty model-implied pairwise cells.",
         call. = FALSE)
  }
  root <- sqrt(1 - rho^2)
  # Mass of the conditioning interval on the other variable, at an edge.
  strip <- function(edge, other_lo, other_hi) {
    out <- numeric(length(edge))
    ok <- is.finite(edge)
    if (!any(ok)) return(out)
    z <- function(t) {
      w <- (t[ok] - rho * edge[ok]) / root
      w[t[ok] == Inf] <- Inf
      w[t[ok] == -Inf] <- -Inf
      stats::pnorm(w)
    }
    out[ok] <- stats::dnorm(edge[ok]) * (z(other_hi) - z(other_lo))
    out
  }
  side <- function(code, tau, edge_lo, edge_hi, other_lo, other_hi) {
    n <- length(code)
    out <- matrix(0, n, length(tau))
    upper <- strip(edge_hi, other_lo, other_hi)
    lower <- strip(edge_lo, other_lo, other_hi)
    for (a in seq_along(tau)) {
      out[code == a, a] <- upper[code == a]
      out[code == a + 1L, a] <- -lower[code == a + 1L]
    }
    out / probability
  }
  list(
    rho = bvn_rect(lo_u, hi_u, lo_v, hi_v, rho, derivative = TRUE) / probability,
    tau_u = side(code_u, tau_u, lo_u, hi_u, lo_v, hi_v),
    tau_v = side(code_v, tau_v, lo_v, hi_v, lo_u, hi_u)
  )
}

#' Analytic per-case influence of the ordinal first stage.
#'
#' Returns the `n` by `m` influence function of the moment vector
#' `[thresholds; lower-triangle polychorics]`, groupwise centered, whose
#' empirical covariance is the first-stage Gamma.
#' @keywords internal
ordinal_moment_influence <- function(X, thresholds, rho, parts = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  sizes <- vapply(thresholds, length, integer(1))
  offsets <- c(0L, cumsum(sizes))
  n_threshold <- sum(sizes)
  pairs <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  n_pair <- nrow(pairs)

  TH <- matrix(0, n, n_threshold)
  for (j in seq_len(p)) {
    TH[, offsets[j] + seq_len(sizes[j])] <-
      threshold_influence(X[, j], thresholds[[j]])
  }

  # Each correlation is estimated at the already-estimated thresholds, so its
  # influence carries theirs through the cross term:
  #
  #   phi_rho = (psi_rho - A21 %*% phi_tau) / A22
  #
  # with A22 and A21 the pair likelihood's own information quantities, taken as
  # outer products, which the information identity licenses within the pair.
  COR <- matrix(0, n, n_pair)
  detail <- if (parts) vector("list", n_pair) else NULL
  for (k in seq_len(n_pair)) {
    u <- pairs[k, 1L]
    v <- pairs[k, 2L]
    ps <- pair_scores(X[, u], X[, v], thresholds[[u]], thresholds[[v]], rho[k])
    a22 <- mean(ps$rho^2)
    c_u <- crossprod(ps$tau_u, ps$rho) / n
    c_v <- crossprod(ps$tau_v, ps$rho) / n
    carry <- TH[, offsets[u] + seq_len(sizes[u]), drop = FALSE] %*% c_u +
      TH[, offsets[v] + seq_len(sizes[v]), drop = FALSE] %*% c_v
    COR[, k] <- (ps$rho - drop(carry)) / a22
    if (parts) {
      detail[[k]] <- list(u = u, v = v, score = ps$rho, tau_u = ps$tau_u,
                          tau_v = ps$tau_v, a22 = a22, c_u = c_u, c_v = c_v)
    }
  }

  influence <- cbind(TH, COR)
  influence <- sweep(influence, 2L, colMeans(influence), "-")
  if (!parts) return(influence)
  list(influence = influence, TH = TH, detail = detail, pairs = pairs,
       sizes = sizes, offsets = offsets, n_threshold = n_threshold)
}

# Derivative of the first-stage weight matrix along the empirical measure.
#
# Perturbing the empirical measure toward a case moves Gamma = E[phi phi'] three
# ways, and the third is easy to miss:
#
#   1. the outer product itself,            phi_i phi_i' - Gamma
#   2. the estimated thresholds and
#      correlations,                        dphi/dkappa . phi_i
#   3. the bread inside phi                 a22 and the carry coefficients are
#                                           themselves sample means, so the
#                                           measure moves them directly
#
# Channel 3 is why a kappa-only derivative sits at a stable wrong ratio rather
# than converging. All three are below; only channel 2 needs finite differences,
# and those are taken on a closed-form function of kappa, not a re-solved
# polychoric optimizer.
ordinal_weight_influence <- function(X, thresholds, rho, step = 1e-5) {
  n <- nrow(X)
  p <- length(thresholds)
  built <- ordinal_moment_influence(X, thresholds, rho, parts = TRUE)
  phi <- built$influence
  gamma <- crossprod(phi) / n
  m <- ncol(phi)
  offsets <- built$offsets
  sizes <- built$sizes
  n_threshold <- built$n_threshold

  # Channels 1 and 3, in closed form. E[phi_k dphi_k/deps_i] for the threshold
  # block is zero: a threshold influence has no bread and no sample mean in it.
  cross <- matrix(0, n, m)
  for (k in seq_along(built$detail)) {
    d <- built$detail[[k]]
    column <- n_threshold + k
    phi_k <- phi[, column]
    iu <- offsets[d$u] + seq_len(sizes[d$u])
    iv <- offsets[d$v] + seq_len(sizes[d$v])
    # E_j[phi_tau,j phi_rho,j], the weights the moved carry coefficients hit.
    g_u <- crossprod(built$TH[, iu, drop = FALSE], phi_k) / n
    g_v <- crossprod(built$TH[, iv, drop = FALSE], phi_k) / n
    # Von Mises derivatives of the sample means that make up the bread.
    d_cu <- sweep(d$tau_u * d$score, 2L, drop(d$c_u), "-")
    d_cv <- sweep(d$tau_v * d$score, 2L, drop(d$c_v), "-")
    d_a22 <- d$score^2 - d$a22
    cross[, column] <- 2 * (
      -(drop(d_cu %*% g_u) + drop(d_cv %*% g_v)) / d$a22 -
        gamma[column, column] * d_a22 / d$a22)
  }

  # Channel 2. kappa moves by phi_i, so accumulate dGamma/dkappa_c weighted by
  # each case's own influence in that coordinate.
  kappa <- c(unlist(thresholds, use.names = FALSE), rho)
  unpack <- function(k) {
    list(thresholds = lapply(seq_len(p),
                             function(j) k[offsets[j] + seq_len(sizes[j])]),
         rho = k[-seq_len(n_threshold)])
  }
  kappa_channel <- matrix(0, m, m)
  quadratic <- function(k) {
    u <- unpack(k)
    colSums(ordinal_moment_influence(X, u$thresholds, u$rho)^2) / n
  }
  for (cc in seq_len(m)) {
    h <- step * max(1, abs(kappa[cc]))
    shift <- function(delta) quadratic(`[<-`(kappa, cc, kappa[cc] + delta))
    # Richardson extrapolation. The integrand is closed form and smooth, so the
    # O(h^4) estimate is essentially step-free, which matters because the step
    # is user-visible and results should not move with it.
    coarse <- (shift(h) - shift(-h)) / (2 * h)
    fine <- (shift(h / 2) - shift(-h / 2)) / h
    kappa_channel[, cc] <- (4 * fine - coarse) / 3
  }

  diagonal <- sweep(phi^2, 2L, diag(gamma), "-") + cross + phi %*% t(kappa_channel)
  list(gamma = gamma, diagonal = diagonal, phi = phi)
}
