#!/usr/bin/env Rscript

# Seeded, non-CRAN validation for the two FIML information conventions.
# Run from the package root. Environment variables:
#   NREP=500 N=400 NCORES=8

suppressPackageStartupMessages({
  library(lavaan)
  library(semTests)
  library(parallel)
})

nrep <- as.integer(Sys.getenv("NREP", "500"))
n <- as.integer(Sys.getenv("N", "400"))
ncores <- as.integer(Sys.getenv(
  "NCORES",
  as.character(min(8L, max(1L, parallel::detectCores() - 2L)))
))

population <- "
  f1 =~ .7*y1 + .7*y2 + .7*y3
  f2 =~ .7*y4 + .7*y5 + .7*y6
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f1 ~~ .4*f2
  y1 ~~ .51*y1
  y2 ~~ .51*y2
  y3 ~~ .51*y3
  y4 ~~ .51*y4
  y5 ~~ .51*y5
  y6 ~~ .51*y6
"

h1_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
"

h0_model <- "
  f1 =~ a*y1 + a*y2 + a*y3
  f2 =~ b*y4 + b*y5 + b*y6
"

make_data <- function(scenario, iteration) {
  set.seed(202607170L + 10000L * match(
    scenario, c("complete_normal", "mcar_normal", "mar_normal", "mar_t5")
  ) + iteration)
  data <- simulateData(
    population, sample.nobs = n, model.type = "cfa", std.lv = TRUE
  )
  if (scenario == "mar_t5") {
    radial <- sqrt(3 / stats::rchisq(n, df = 5))
    data[] <- lapply(data, function(x) x * radial)
  }
  if (scenario == "mcar_normal") {
    miss1 <- stats::runif(n) < .30
    miss2 <- stats::runif(n) < .30
    data[miss1, c("y1", "y2", "y3")] <- NA_real_
    data[miss2, c("y4", "y5")] <- NA_real_
  } else if (scenario %in% c("mar_normal", "mar_t5")) {
    z <- as.numeric(scale(data$y6))
    probability1 <- stats::plogis(stats::qlogis(.30) + .65 * z)
    probability2 <- stats::plogis(stats::qlogis(.30) - .65 * z)
    miss1 <- stats::runif(n) < probability1
    miss2 <- stats::runif(n) < probability2
    data[miss1, c("y1", "y2", "y3")] <- NA_real_
    data[miss2, c("y4", "y5")] <- NA_real_
  }
  data
}

one <- function(scenario, iteration) {
  data <- make_data(scenario, iteration)
  out <- rep(NA_real_, 16L)
  names(out) <- as.vector(outer(
    c("observed_delta", "observed_exact", "lavaan_delta", "lavaan_exact"),
    c("sb", "ss", "all", "peba4"),
    paste, sep = "_"
  ))
  try({
    m1 <- cfa(
      h1_model, data, missing = "fiml", estimator = "MLR", std.lv = TRUE
    )
    m0 <- cfa(
      h0_model, data, missing = "fiml", estimator = "MLR", std.lv = TRUE
    )
    if (!lavInspect(m0, "converged") || !lavInspect(m1, "converged")) {
      return(out)
    }
    df <- as.integer(fitMeasures(m0, "df") - fitMeasures(m1, "df"))
    statistic <- as.numeric(
      fitMeasures(m0, "chisq") - fitMeasures(m1, "chisq")
    )
    for (convention in c("observed", "lavaan")) {
      for (A.method in c("delta", "exact")) {
        lambda <- semTests:::fiml_lambdas_nested(
          m0, m1, df,
          A.method = A.method,
          fiml.convention = convention
        )$ug_biased
        prefix <- paste(convention, A.method, sep = "_")
        out[paste0(prefix, "_sb")] <-
          semTests:::trad_pvalue(df, statistic, lambda, "sb")
        out[paste0(prefix, "_ss")] <-
          semTests:::trad_pvalue(df, statistic, lambda, "ss")
        out[paste0(prefix, "_all")] <-
          semTests:::trad_pvalue(df, statistic, lambda, "all")
        out[paste0(prefix, "_peba4")] <-
          semTests:::peba_pvalue(statistic, lambda, 4L)
      }
    }
  }, silent = TRUE)
  out
}

scenarios <- c("complete_normal", "mcar_normal", "mar_normal", "mar_t5")
results <- lapply(scenarios, function(scenario) {
  values <- parallel::mclapply(
    seq_len(nrep), function(i) one(scenario, i),
    mc.cores = ncores, mc.preschedule = TRUE
  )
  values <- do.call(rbind, values)
  valid <- colSums(is.finite(values))
  rejection <- vapply(c(.01, .05, .10), function(alpha) {
    colMeans(values < alpha, na.rm = TRUE)
  }, numeric(ncol(values)))
  data.frame(
    scenario = scenario,
    method = colnames(values),
    valid = valid,
    alpha_01 = rejection[, 1L],
    alpha_05 = rejection[, 2L],
    alpha_10 = rejection[, 3L],
    row.names = NULL
  )
})
results <- do.call(rbind, results)

dir.create("workspace", showWarnings = FALSE)
utils::write.csv(
  results, "workspace/fiml-validation-simulation-results.csv",
  row.names = FALSE
)
print(results, row.names = FALSE)
