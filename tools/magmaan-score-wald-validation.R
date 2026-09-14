#!/usr/bin/env Rscript
# Optional planning-stage primitive checks; not a certification of a lavaan adapter.
# Run from the semTests root. magmaan is deliberately not a package dependency.
stopifnot(requireNamespace("magmaan", quietly = TRUE),
          requireNamespace("pkgload", quietly = TRUE))
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
rows <- list()
check <- function(label, actual, expected, tolerance = 2e-6) {
  stopifnot(length(actual) == length(expected), length(actual) > 0L,
            all(is.finite(actual)), all(is.finite(expected)))
  gap <- max(abs(actual - expected))
  rows[[length(rows) + 1L]] <<- data.frame(check = label, gap, tolerance)
  if (gap > tolerance) stop(label, ": gap ", gap)
}
model <- "visual =~ x1+x2+x3
          textual =~ x4+x5+x6
          speed =~ x7+x8+x9"
original <- lavaan::HolzingerSwineford1939
for (estimator in c("ML", "FIML")) {
  dat <- original
  if (estimator == "FIML") {
    set.seed(20260912)
    dat[sample.int(nrow(dat), 60), "x2"] <- NA
    dat[sample.int(nrow(dat), 45), "x5"] <- NA
  }
  fit <- magmaan::magmaan(model, dat, estimator = estimator)
  stopifnot(isTRUE(fit$converged))
  for (sensitivity in c("expected", "observed")) {
    z <- magmaan::global_score_flip_test(
      fit, if (estimator == "FIML") NULL else dat,
      n_flips = 19L, seed = 1, sensitivity = sensitivity)
    # Few flips only exercise the wrapper; no multiplier calibration claim.
    s <- z$projected_score
    H <- z$projected_metric
    J <- z$projected_meat
    invroot <- eigen(H, symmetric = TRUE)
    stopifnot(all(invroot$values > 0))
    W <- invroot$vectors %*% diag(1 / sqrt(invroot$values)) %*%
      t(invroot$vectors)
    lambda <- eigen(W %*% J %*% W, symmetric = TRUE)$values
    prefix <- paste(estimator, sensitivity)
    check(paste(prefix, "score quadratic"), drop(crossprod(s, solve(H, s))),
          z$statistic_effective)
    check(paste(prefix, "spectrum"), sort(lambda), sort(z$eigenvalues))
    check(paste(prefix, "sandwich quadratic"), drop(crossprod(s, solve(J, s))),
          z$statistic_sandwich)
    check(paste(prefix, "mixture p"),
          semTests:::pvalue_all(z$statistic_effective, lambda), z$p_mixture)
    check(paste(prefix, "SB p"),
          pchisq(z$statistic_effective / mean(lambda), z$df, lower.tail = FALSE),
          z$p_mean_scaled)
    result <- semTests::score_test_matrices(s, H, J)
    studentized <- semTests::score_test_matrices(s, score_covariance = J)
    check(paste(prefix, "public score"), result$statistic, z$statistic_effective)
    check(paste(prefix, "public spectrum"), sort(result$eigenvalues), sort(z$eigenvalues))
    check(paste(prefix, "public mixture"), as.numeric(semTests::pvalues(result)), z$p_mixture)
    check(paste(prefix, "public SB"), as.numeric(semTests::pvalues(result, "SB")), z$p_mean_scaled)
    check(paste(prefix, "public studentized"), studentized$statistic, z$statistic_sandwich)
    check(paste(prefix, "public studentized p"), studentized$p_value, z$p_sandwich)
  }
}
theta <- c(0.4, -0.2, 0.7)
V <- matrix(c(2, .2, .1, .2, 1, .3, .1, .3, 1.5), 3) / 200
R <- rbind(c(1, -1, 0), c(0, 1, -1))
q <- c(.1, -.3)
z <- magmaan::magmaan_core$infer_wald_test_theta(theta, R, V, q)
contrast <- drop(R %*% theta) - q
reference <- drop(crossprod(contrast, solve(R %*% V %*% t(R), contrast)))
check("linear Wald statistic", z$chi2, reference)
check("linear Wald p", z$pvalue, pchisq(reference, nrow(R), lower.tail = FALSE))
wald <- semTests::wald_test_matrices(theta, R, V, q)
check("public Wald statistic", wald$statistic, z$chi2)
check("public Wald p", as.numeric(semTests::pvalues(wald)), z$pvalue)
dir.create("workspace", showWarnings = FALSE)
write.csv(do.call(rbind, rows), "workspace/magmaan-score-wald-planning.csv",
          row.names = FALSE)
print(do.call(rbind, rows), row.names = FALSE)
cat(length(rows), "primitive checks passed; adapter and calibration remain untested.\n")
