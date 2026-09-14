#!/usr/bin/env Rscript
# Fixed fitted objects: compare spectrum algorithms without timing estimation.
# REPS=3 P=9,18,30 Rscript tools/residual-spectrum-timing.R
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
reps <- as.integer(Sys.getenv("REPS", "3"))
sizes <- as.integer(strsplit(Sys.getenv("P", "9,18,30"), ",", fixed = TRUE)[[1]])
stopifnot(!is.na(reps), reps > 0, all(is.finite(sizes)),
          all(sizes >= 6 & sizes %% 3 == 0))
elapsed <- function(f) {
  f()
  unname(system.time(for (i in seq_len(reps)) f())[["elapsed"]]) * 1000 / reps
}
rows <- lapply(sizes, function(p) {
  set.seed(9212 + p)
  factors <- matrix(rnorm(500 * 3), 500)
  X <- .7 * factors[, rep(1:3, each = p / 3)] + matrix(rnorm(500 * p), 500)
  X <- as.data.frame(X)
  names(X) <- paste0("x", seq_len(p))
  model <- paste(vapply(1:3, function(k) {
    paste0("f", k, " =~ ", paste(names(X)[((k - 1) * p / 3 + 1):(k * p / 3)],
                                  collapse = " + "))
  }, character(1)), collapse = "\n")
  fit <- lavaan::cfa(model, X, estimator = "MLM")
  stopifnot(isTRUE(lavaan::lavInspect(fit, "converged")))
  df <- semTests:::fit_df(fit)
  reference <- function() semTests:::lavaan_lambdas(fit, df)
  reduced <- function() semTests:::single_model_reduced_lambdas(fit, df)
  gap <- max(abs(reference()$ug_biased - reduced()$ug_biased))
  stopifnot(is.finite(gap), gap < 1e-8)
  data.frame(p, moments = p * (p + 1) / 2, df = as.integer(df), gap,
             full_ms = elapsed(reference), reduced_ms = elapsed(reduced))
})
out <- do.call(rbind, rows)
out$speedup <- out$full_ms / out$reduced_ms
dir.create("workspace", showWarnings = FALSE)
write.csv(out, "workspace/residual-spectrum-timing.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "workspace/residual-spectrum-session.txt")
print(out, row.names = FALSE)
