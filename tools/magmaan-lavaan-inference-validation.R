#!/usr/bin/env Rscript
# Independent fits and independently extracted ML global-score ingredients.
# The Wald checks intentionally share lavaan's fitted covariance with the C++
# primitive; they do not certify a new covariance estimator.
stopifnot(requireNamespace("magmaan", quietly = TRUE),
          requireNamespace("pkgload", quietly = TRUE))
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
model <- "visual =~ x1+x2+x3
          textual =~ x4+x5+x6
          speed =~ x7+x8+x9"
dat <- lavaan::HolzingerSwineford1939
shifted <- dat
shifted$x1 <- shifted$x1 - mean(shifted$x1) + .2
cases <- list(
  covariance = list(model = model, data = dat, meanstructure = FALSE),
  means = list(model = model, data = dat, meanstructure = TRUE),
  fixed_mean = list(model = paste(model, "x1 ~ 0*1", sep = "\n"),
                    data = shifted, meanstructure = TRUE),
  unequal_groups = list(model = model, data = dat, meanstructure = TRUE,
                        group = "school", group.equal = "loadings"),
  random_x = list(model = "f =~ x1+x2+x3\nf ~ ageyr", data = dat,
                  meanstructure = TRUE)
)
rows <- list()
check <- function(case, quantity, actual, expected, tolerance) {
  stopifnot(length(actual) == length(expected), length(actual) > 0L,
            all(is.finite(actual)), all(is.finite(expected)))
  gap <- max(abs(actual - expected))
  rows[[length(rows) + 1L]] <<- data.frame(case, quantity, gap, tolerance)
  if (gap > tolerance) stop(case, " ", quantity, ": gap ", gap)
}
for (case in names(cases)) {
  args <- cases[[case]]
  fit <- do.call(lavaan::cfa, c(args, list(estimator = "MLM", fixed.x = FALSE)))
  spec <- magmaan::model_spec(args$model, meanstructure = args$meanstructure,
    fixed_x = FALSE, group = args$group, group_equal = args$group.equal,
    group_labels = if (!is.null(args$group)) unique(as.character(args$data[[args$group]])) else NULL)
  mag <- magmaan::magmaan(spec, args$data, estimator = "ML")
  stopifnot(isTRUE(lavaan::lavInspect(fit, "converged")), isTRUE(mag$converged))
  reference <- magmaan::global_score_flip_test(mag, args$data, n_flips = 1L)
  score <- semTests::score_test(fit)
  sandwich <- semTests::score_test(fit, metric = "sandwich")
  # Separate optimizers: tolerances concern endpoint differences, not just
  # floating-point algebra. One flip only enables the deterministic comparators.
  check(case, "df", score$df, reference$df, 0)
  check(case, "score", score$statistic, reference$statistic_effective, 2e-3)
  check(case, "spectrum", sort(score$eigenvalues), sort(reference$eigenvalues), 5e-4)
  check(case, "mixture", score$p_value, reference$p_mixture, 1e-4)
  check(case, "SB", as.numeric(semTests::pvalues(score, "SB")), reference$p_mean_scaled, 1e-4)
  check(case, "sandwich", sandwich$statistic, reference$statistic_sandwich, 2e-3)
  check(case, "sandwich p", sandwich$p_value, reference$p_sandwich, 1e-4)
  theta <- lavaan::coef(fit)
  R <- matrix(0, 1, length(theta)); R[1, 1:2] <- c(1, -1)
  wald <- semTests::wald_test(fit, R)
  native <- magmaan::magmaan_core$infer_wald_test_theta(
    as.numeric(theta), R, as.matrix(lavaan::lavInspect(fit, "vcov")), 0)
  check(case, "Wald supplied-covariance statistic", wald$statistic, native$chi2, 1e-8)
  check(case, "Wald supplied-covariance p", wald$p_value, native$pvalue, 1e-8)
}
out <- do.call(rbind, rows)
dir.create("workspace", showWarnings = FALSE)
write.csv(out, "workspace/magmaan-lavaan-inference.csv", row.names = FALSE)
files <- list.files(system.file(package = "magmaan"), pattern = "\\.(rdb|so)$",
                    recursive = TRUE, full.names = TRUE)
capture.output(sessionInfo(), tools::md5sum(files),
               file = "workspace/magmaan-lavaan-inference-session.txt")
print(out, row.names = FALSE)
cat(nrow(out), "adapter checks passed; no simulation calibration claim.\n")
