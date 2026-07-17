#!/usr/bin/env Rscript

# Seeded, non-CRAN validation for categorical single and nested tests.
# Run from the package root. Environment variables:
#   CORE_REPS=500 STRESS_REPS=300 N=500 NCORES=8

suppressPackageStartupMessages({
  library(lavaan)
  library(semTests)
  library(parallel)
})

core_reps <- as.integer(Sys.getenv("CORE_REPS", "500"))
stress_reps <- as.integer(Sys.getenv("STRESS_REPS", "300"))
n <- as.integer(Sys.getenv("N", "500"))
ncores <- as.integer(Sys.getenv(
  "NCORES",
  as.character(min(8L, max(1L, parallel::detectCores() - 2L)))
))

population <- "
  f1 =~ .7*y1 + .7*y2 + .7*y3
  f2 =~ .7*y4 + .7*y5 + .7*y6
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f1 ~~ .35*f2
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

h0_multigroup <- "
  f1 =~ c(a1, a2)*y1 + c(a1, a2)*y2 + c(a1, a2)*y3
  f2 =~ c(b1, b2)*y4 + c(b1, b2)*y5 + c(b1, b2)*y6
"

scenario_table <- data.frame(
  scenario = c(
    "wlsmv_complete", "wlsmv_pairwise_mcar",
    "ulsmv_complete", "ulsmv_pairwise_mcar",
    "wlsmv_binary_skew", "wlsmv_mixed",
    "wlsmv_multigroup", "wlsmv_pairwise_mar"
  ),
  estimator = c(
    "WLSMV", "WLSMV", "ULSMV", "ULSMV",
    "WLSMV", "WLSMV", "WLSMV", "WLSMV"
  ),
  reps = c(
    rep(core_reps, 4L), rep(stress_reps, 4L)
  ),
  stringsAsFactors = FALSE
)

make_data <- function(scenario, iteration) {
  scenario_index <- match(scenario, scenario_table$scenario)
  set.seed(202607190L + scenario_index * 10000L + iteration)
  data <- simulateData(
    population, sample.nobs = n, model.type = "cfa", std.lv = TRUE
  )
  mar_driver <- as.numeric(scale(data$y6))

  ordered_variables <- if (scenario == "wlsmv_mixed") {
    paste0("y", 1:3)
  } else {
    paste0("y", 1:6)
  }
  cutpoints <- if (scenario == "wlsmv_binary_skew") {
    c(-Inf, stats::qnorm(.80), Inf)
  } else {
    c(-Inf, -.5, .5, Inf)
  }
  for (variable in ordered_variables) {
    data[[variable]] <- ordered(cut(
      data[[variable]], breaks = cutpoints,
      include.lowest = TRUE
    ))
  }

  if (grepl("pairwise_mcar$", scenario)) {
    for (variable in paste0("y", 1:6)) {
      data[stats::runif(n) < .25, variable] <- NA
    }
  } else if (scenario == "wlsmv_pairwise_mar") {
    for (j in seq_len(5L)) {
      probability <- stats::plogis(
        stats::qlogis(.25) + (-1)^j * .7 * mar_driver
      )
      data[stats::runif(n) < probability, paste0("y", j)] <- NA
    }
  }
  if (scenario == "wlsmv_multigroup") {
    data$group <- rep(c("g1", "g2"), length.out = n)
  }

  list(
    data = data,
    ordered = ordered_variables,
    missing = if (grepl("pairwise", scenario)) "pairwise" else "listwise",
    group = if (scenario == "wlsmv_multigroup") "group" else NULL
  )
}

one <- function(scenario, estimator, iteration) {
  fixture <- make_data(scenario, iteration)
  output <- rep(NA_real_, 8L)
  names(output) <- c(
    paste0("single_", c("sb", "ss", "all", "peba4")),
    paste0("nested_", c("sb", "ss", "all", "peba4"))
  )
  try({
    common <- list(
      data = fixture$data,
      ordered = fixture$ordered,
      estimator = estimator,
      std.lv = TRUE,
      missing = fixture$missing,
      group = fixture$group
    )
    m1 <- do.call(lavaan::cfa, c(list(h1_model), common))
    m0 <- do.call(
      lavaan::cfa,
      c(list(
        if (scenario == "wlsmv_multigroup") h0_multigroup else h0_model
      ), common)
    )
    if (!lavInspect(m0, "converged") || !lavInspect(m1, "converged")) {
      return(output)
    }
    output[seq_len(4L)] <- unname(pvalues(
      m1, c("SB", "SS", "ALL", "PEBA4")
    ))
    output[4L + seq_len(4L)] <- unname(pvalues_nested(
      m0, m1, tests = c("SB", "SS", "ALL", "PEBA4"),
      method = "2000", A.method = "delta"
    ))
  }, silent = TRUE)
  output
}

results <- lapply(seq_len(nrow(scenario_table)), function(row) {
  specification <- scenario_table[row, ]
  values <- parallel::mclapply(
    seq_len(specification$reps),
    function(iteration) {
      one(specification$scenario, specification$estimator, iteration)
    },
    mc.cores = ncores,
    mc.preschedule = TRUE
  )
  values <- do.call(rbind, values)
  valid <- colSums(is.finite(values))
  rejection <- vapply(c(.01, .05, .10), function(alpha) {
    colMeans(values < alpha, na.rm = TRUE)
  }, numeric(ncol(values)))
  data.frame(
    scenario = specification$scenario,
    method = colnames(values),
    requested_reps = specification$reps,
    valid = valid,
    sample_n = n,
    alpha_01 = rejection[, 1L],
    alpha_05 = rejection[, 2L],
    alpha_10 = rejection[, 3L],
    lavaan_version = as.character(utils::packageVersion("lavaan")),
    semtests_version = as.character(utils::packageVersion("semTests")),
    row.names = NULL
  )
})
results <- do.call(rbind, results)

dir.create("workspace", showWarnings = FALSE)
utils::write.csv(
  results, "workspace/categorical-validation-simulation-results.csv",
  row.names = FALSE
)
print(results, row.names = FALSE)
