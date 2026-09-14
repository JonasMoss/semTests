#!/usr/bin/env Rscript
# Compare covariance projection only, excluding fitting and eigenvalues.
# OPENBLAS_NUM_THREADS=1 REPS=5 Rscript tools/projected-gamma-timing.R
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
reps <- as.integer(Sys.getenv("REPS", "5"))
sizes <- as.integer(strsplit(Sys.getenv("P", "9,18,30"), ",", fixed = TRUE)[[1]])
stopifnot(!is.na(reps), reps > 0L, all(is.finite(sizes)),
          all(sizes >= 6 & sizes %% 3 == 0))
elapsed <- function(f) {
  f()
  unname(system.time(for (i in seq_len(reps)) f())[["elapsed"]]) * 1000 / reps
}
rows <- list()
for (p in sizes) {
  set.seed(991 + p)
  latent <- matrix(rnorm(500 * 3), 500)
  X <- as.data.frame(.7 * latent[, rep(1:3, each = p / 3)] +
                      matrix(rnorm(500 * p), 500))
  names(X) <- paste0("x", seq_len(p))
  model <- paste(vapply(1:3, function(k) {
    paste0("f", k, " =~ ", paste(names(X)[((k - 1) * p / 3 + 1):(k * p / 3)],
                                  collapse = " + "))
  }, character(1)), collapse = "\n")
  fit <- lavaan::cfa(model, X, estimator = "ML", se = "none", test = "standard")
  stopifnot(isTRUE(lavaan::lavInspect(fit, "converged")))
  cached <- fit
  cached@SampleStats@NACOV <- semTests:::gamma_from_lavaan(fit)
  B <- semTests:::residual_metric_factor(lavaan::lavTech(fit, "WLS.V"),
                                        semTests:::model_effective_delta(fit),
                                        semTests:::fit_df(fit))
  for (r in c(3L, ncol(B))) {
    F <- B[, seq_len(r), drop = FALSE]
    full <- function() semTests:::projected_gamma(fit, F, source = "materialized")
    raw <- function() semTests:::projected_gamma(fit, F, source = "raw")
    reuse <- function() semTests:::projected_gamma(cached, F)
    gap <- max(abs(full() - raw()))
    stopifnot(is.finite(gap), gap < 1e-8)
    rows[[length(rows) + 1L]] <- data.frame(
      p, n = nrow(X), moments = nrow(F), projected = r, gap,
      full_uncached_ms = elapsed(full), tiled_raw_ms = elapsed(raw),
      cached_ms = elapsed(reuse))
  }
}
out <- do.call(rbind, rows)
out$raw_speedup_over_uncached <- out$full_uncached_ms / out$tiled_raw_ms
dir.create("workspace", showWarnings = FALSE)
write.csv(out, "workspace/projected-gamma-timing.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "workspace/projected-gamma-session.txt")
print(out, row.names = FALSE)
