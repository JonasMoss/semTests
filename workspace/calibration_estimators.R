#!/usr/bin/env Rscript
# Calibration harness for the broadened estimator support (the Phase-0 gate).
#
# Type-I rejection rates at nominal alpha for the eigenvalue p-values, under a
# CORRECTLY specified 2-factor CFA, across estimators / data types. A final
# misspecified cell checks that power is non-degenerate. This is a workspace
# script (build-ignored); it is NOT a testthat test.
#
#   Rscript workspace/calibration_estimators.R [reps] [n1,n2,...]
#
# Defaults: reps = 1000, n = 300,1000. Monte-Carlo SE ~ sqrt(a(1-a)/reps);
# at reps=1000, alpha=.05 -> SE ~ .007. Well-calibrated = empirical within a
# few SE of nominal.

suppressMessages({
  library(lavaan)
  pkgload::load_all(getwd(), quiet = TRUE)
})

args  <- commandArgs(trailingOnly = TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
ns    <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(300L, 1000L)
alphas <- c(0.01, 0.05, 0.10)
tests <- c("PEBA4", "PALL", "SS", "SB")        # eigenvalue corrections + SB reference

# --- population: 2 factors, 3 indicators each, standardised loadings 0.7 ------
loads  <- 0.7
resid  <- 1 - loads^2
pop <- paste(
  sprintf("F1 =~ %g*x1 + %g*x2 + %g*x3", loads, loads, loads),
  sprintf("F2 =~ %g*x4 + %g*x5 + %g*x6", loads, loads, loads),
  "F1 ~~ 1*F1", "F2 ~~ 1*F2", "F1 ~~ 0.3*F2",
  paste(sprintf("x%d ~~ %g*x%d", 1:6, resid, 1:6), collapse = "\n"),
  sep = "\n"
)
fitmod  <- "F1 =~ x1 + x2 + x3\nF2 =~ x4 + x5 + x6"   # correctly specified
missmod <- "F =~ x1 + x2 + x3 + x4 + x5 + x6"          # 1-factor: misspecified

fit_one <- function(estimator, dat, model) {
  switch(estimator,
    ML    = cfa(model, dat, estimator = "MLM"),
    GLS   = cfa(model, dat, estimator = "GLS"),
    ULS   = cfa(model, dat, estimator = "ULS", test = "satorra.bentler"),
    MLR   = cfa(model, dat, estimator = "MLR"),
    FIML  = {
      d <- dat; d[matrix(stats::runif(nrow(dat) * 6) < 0.15, nrow(dat), 6)] <- NA
      cfa(model, d, missing = "fiml", estimator = "MLR")
    },
    WLSMV = {
      d <- as.data.frame(lapply(dat, function(z)
        ordered(cut(z, breaks = c(-Inf, -0.6, 0.6, Inf)))))
      cfa(model, d, ordered = names(d))
    }
  )
}

estimators <- c("ML", "GLS", "ULS", "MLR", "FIML", "WLSMV")

run_cell <- function(estimator, n, model, seed0) {
  rej <- matrix(0, length(tests), length(alphas),
                dimnames = list(tests, paste0("a", alphas)))
  valid <- 0L
  for (r in seq_len(reps)) {
    set.seed(seed0 + r)
    p <- tryCatch(suppressWarnings({
      dat <- simulateData(pop, sample.nobs = n)
      fit <- fit_one(estimator, dat, model)
      unlist(lapply(tests, function(t) as.numeric(pvalues(fit, t))))
    }), error = function(e) NULL)
    if (is.null(p) || length(p) != length(tests) || anyNA(p)) next
    valid <- valid + 1L
    for (j in seq_along(alphas)) rej[, j] <- rej[, j] + (p < alphas[j])
  }
  data.frame(
    estimator = estimator, n = n,
    test = rep(tests, length(alphas)),
    alpha = rep(alphas, each = length(tests)),
    reject = as.vector(rej) / max(valid, 1L),
    valid = valid,
    row.names = NULL
  )
}

cat(sprintf("Calibration: reps=%d, n=%s, tests=%s\n\n",
            reps, paste(ns, collapse = ","), paste(tests, collapse = ",")))

## --- Type-I error (correct model) --------------------------------------------
out <- list()
for (n in ns) for (est in estimators) {
  cat(sprintf("[type-I] %-6s n=%-5d ...\n", est, n))
  out[[length(out) + 1L]] <- run_cell(est, n, fitmod, seed0 = 10000L * which(ns == n))
}
typeI <- do.call(rbind, out)

cat("\n===== Type-I rejection rates (correct model; target = alpha) =====\n")
print(reshape(typeI[c("estimator", "n", "test", "alpha", "reject")],
              idvar = c("estimator", "n", "test"), timevar = "alpha",
              direction = "wide"), row.names = FALSE, digits = 3)

## --- Power (misspecified: fit 1-factor) at n = max --------------------------
np <- max(ns)
cat(sprintf("\n===== Power (1-factor fit to 2-factor data; n=%d) =====\n", np))
pw <- do.call(rbind, lapply(estimators, function(est) {
  cat(sprintf("[power]  %-6s ...\n", est))
  run_cell(est, np, missmod, seed0 = 99000L)
}))
pw5 <- pw[pw$alpha == 0.05, c("estimator", "test", "reject", "valid")]
print(pw5, row.names = FALSE, digits = 3)

saveRDS(list(typeI = typeI, power = pw, reps = reps, ns = ns),
        "workspace/calibration_estimators.Rds")
cat("\nSaved workspace/calibration_estimators.Rds\n")
