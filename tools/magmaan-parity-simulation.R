#!/usr/bin/env Rscript

# Seeded, non-CRAN SIMULATION parity check: semTests vs magmaan across many
# randomly generated datasets per estimator cell. This complements the
# single-dataset tools/magmaan-validation.R gate by confirming that the two
# independent implementations still agree once each fits its own draw of the
# data -- i.e. that the single-example parity is not a fixed-dataset accident.
#
# It measures agreement (max |p_semTests - p_magmaan| over the test family per
# draw), NOT calibration / rejection rates. Run from the package root:
#
#   REPS=200 N=500 NCORES=8 CELLS=all \
#     PARITY_SIM_OUTPUT=workspace/magmaan-parity-simulation.csv \
#     Rscript tools/magmaan-parity-simulation.R
#
# CELLS is a comma-separated filter of cell names (default "all").

if (!requireNamespace("magmaan", quietly = TRUE)) {
  stop("Install magmaan's r-package/ before running this optional check.",
       call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to validate the current semTests source tree.",
       call. = FALSE)
}

# Keep math libraries single-threaded so mclapply forks do not oversubscribe.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(lavaan)
  library(magmaan)
  library(parallel)
  # Load the current working tree, not the installed build -- the same choice
  # tools/magmaan-validation.R makes, so both check the code under development.
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
})

reps <- as.integer(Sys.getenv("REPS", "200"))
n <- as.integer(Sys.getenv("N", "500"))
ncores <- as.integer(Sys.getenv(
  "NCORES",
  as.character(min(8L, max(1L, parallel::detectCores() - 2L)))
))
cells_filter <- Sys.getenv("CELLS", "all")
sim_output <- Sys.getenv("PARITY_SIM_OUTPUT", "")

# Build the test family for a given df, sizing pEBA blocks like the gate does
# (min(4, df)); pEBA/EBA cannot use more blocks than the test has df.
tests_for_df <- function(df) {
  df <- as.integer(df)
  if (length(df) != 1L || is.na(df) || df < 1L) return(character())
  blocks <- max(1L, min(4L, df))
  c("SB", "SS", "ALL", paste0("PEBA", blocks))
}

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

# H0 keeps the first indicator as the scaling marker and constrains only the
# remaining loadings equal -- a clean nested restriction under lavaan's default
# marker identification (matches the gate's HolzingerSwineford H0 pattern).
h0_model <- "
  f1 =~ y1 + a*y2 + a*y3
  f2 =~ y4 + b*y5 + b*y6
"

ordered_names <- paste0("y", seq_len(6L))

continuous_draw <- function(seed) {
  set.seed(seed)
  lavaan::simulateData(
    population, sample.nobs = n, model.type = "cfa", std.lv = TRUE
  )
}

ordinalize <- function(data, vars = ordered_names) {
  cutpoints <- c(-Inf, -.5, .5, Inf)
  for (variable in vars) {
    data[[variable]] <- ordered(cut(
      data[[variable]], breaks = cutpoints, include.lowest = TRUE
    ))
  }
  data
}

nested_df_of <- function(lav_h0, lav_h1) {
  as.integer(
    lavaan::fitmeasures(lav_h0, "df") - lavaan::fitmeasures(lav_h1, "df")
  )
}

# max |p_semTests - p_magmaan| with alignment guards.
p_diff <- function(semtests, magmaan) {
  semtests <- as.numeric(semtests)
  magmaan <- as.numeric(magmaan)
  if (length(semtests) != length(magmaan) || !length(semtests) ||
      any(!is.finite(semtests)) || any(!is.finite(magmaan))) {
    return(NA_real_)
  }
  max(abs(semtests - magmaan))
}

na_pair <- c(single = NA_real_, nested = NA_real_)

# ---- cell runners: each returns c(single = <max diff>, nested = <max diff>) ---

run_continuous_ls <- function(seed, estimator, gamma) {
  out <- na_pair
  try({
    data <- continuous_draw(seed)
    lav_args_h1 <- list(model = h1_model, data = data, estimator = estimator)
    lav_args_h0 <- list(model = h0_model, data = data, estimator = estimator)
    if (identical(estimator, "ULS")) {
      lav_args_h1$test <- lav_args_h0$test <- "satorra.bentler"
    }
    lav_h1 <- do.call(lavaan::cfa, lav_args_h1)
    lav_h0 <- do.call(lavaan::cfa, lav_args_h0)
    if (!lavInspect(lav_h1, "converged") ||
        !lavInspect(lav_h0, "converged")) {
      return(out)
    }
    single_tests <- tests_for_df(lavaan::fitmeasures(lav_h1, "df"))
    nested_tests <- tests_for_df(nested_df_of(lav_h0, lav_h1))
    # Full WLS (ADF) reuses lavaan's H1 weight matrix in magmaan.
    weight <- NULL
    if (identical(estimator, "WLS")) {
      value <- lavaan::lavTech(lav_h1, "WLS.V")
      if (is.list(value)) value <- value[[1L]]
      weight <- as.matrix(value)
    }
    spec_h1 <- magmaan::model_spec(h1_model)
    spec_h0 <- magmaan::model_spec(h0_model)
    mag_h1 <- magmaan::magmaan(
      spec_h1, magmaan::df_to_data(data, spec_h1, scaling = "n-1"),
      estimator = estimator, W = weight
    )
    mag_h0 <- magmaan::magmaan(
      spec_h0, magmaan::df_to_data(data, spec_h0, scaling = "n-1"),
      estimator = estimator, W = weight
    )
    if (!isTRUE(mag_h1$converged) || !isTRUE(mag_h0$converged)) return(out)

    mag_single <- magmaan::fmg_tests(
      mag_h1, tests = single_tests, gamma = gamma, weight = weight
    )
    out["single"] <- p_diff(
      semTests::pvalues(lav_h1, single_tests), mag_single$p_value
    )
    mag_nested <- magmaan::fmg_nested(
      mag_h1, mag_h0, data = data, tests = nested_tests,
      A.method = "delta", gamma = gamma, weight = weight
    )
    out["nested"] <- p_diff(
      semTests::pvalues_nested(
        lav_h0, lav_h1, method = "2000", tests = nested_tests
      ),
      mag_nested$p_value
    )
  }, silent = TRUE)
  out
}

run_fiml <- function(seed) {
  out <- na_pair
  try({
    data <- continuous_draw(seed)
    miss1 <- stats::runif(n) < .25
    miss2 <- stats::runif(n) < .25
    data[miss1, c("y1", "y2", "y3")] <- NA_real_
    data[miss2, c("y4", "y5")] <- NA_real_
    lav_h1 <- lavaan::cfa(
      h1_model, data, missing = "fiml", estimator = "MLR",
      meanstructure = TRUE, std.lv = TRUE
    )
    lav_h0 <- lavaan::cfa(
      h0_model, data, missing = "fiml", estimator = "MLR",
      meanstructure = TRUE, std.lv = TRUE
    )
    if (!lavInspect(lav_h1, "converged") ||
        !lavInspect(lav_h0, "converged")) {
      return(out)
    }
    single_tests <- tests_for_df(lavaan::fitmeasures(lav_h1, "df"))
    nested_tests <- tests_for_df(nested_df_of(lav_h0, lav_h1))
    mag_h1 <- magmaan::magmaan(h1_model, data, estimator = "FIML")
    mag_h0 <- magmaan::magmaan(h0_model, data, estimator = "FIML")
    if (!isTRUE(mag_h1$converged) || !isTRUE(mag_h0$converged)) return(out)

    mag_single <- magmaan::fmg_tests(mag_h1, tests = single_tests)
    out["single"] <- p_diff(
      semTests::pvalues(lav_h1, single_tests, fiml.convention = "observed"),
      mag_single$p_value
    )
    mag_nested <- magmaan::fmg_nested(
      mag_h1, mag_h0, tests = nested_tests, A.method = "delta"
    )
    out["nested"] <- p_diff(
      semTests::pvalues_nested(
        lav_h0, lav_h1, tests = nested_tests,
        A.method = "delta", fiml.convention = "observed"
      ),
      mag_nested$p_value
    )
  }, silent = TRUE)
  out
}

run_categorical <- function(seed, pairwise,
                            lavaan_estimator = "WLSMV",
                            magmaan_estimator = "DWLS") {
  out <- na_pair
  try({
    data <- ordinalize(continuous_draw(seed))
    missing <- if (pairwise) "pairwise" else "listwise"
    if (pairwise) {
      set.seed(seed + 7L)
      for (variable in ordered_names) {
        data[stats::runif(nrow(data)) < .20, variable] <- NA
      }
    }
    lav_h1 <- lavaan::cfa(
      h1_model, data, ordered = ordered_names,
      estimator = lavaan_estimator, missing = missing
    )
    lav_h0 <- lavaan::cfa(
      h0_model, data, ordered = ordered_names,
      estimator = lavaan_estimator, missing = missing
    )
    if (!lavInspect(lav_h1, "converged") ||
        !lavInspect(lav_h0, "converged")) {
      return(out)
    }
    single_tests <- tests_for_df(lavaan::fitmeasures(lav_h1, "df"))
    nested_tests <- tests_for_df(nested_df_of(lav_h0, lav_h1))
    spec_h1 <- magmaan::model_spec(h1_model, ordered = ordered_names)
    spec_h0 <- magmaan::model_spec(h0_model, ordered = ordered_names)
    stats <- magmaan::magmaan_core$data_ordinal_stats_from_df(
      data, spec_h1, ordered = ordered_names, missing = missing,
      pd_gamma = "overlap",
      full_wls_weight = identical(magmaan_estimator, "WLS")
    )
    mag_h1 <- magmaan::magmaan(spec_h1, stats, estimator = magmaan_estimator)
    mag_h0 <- magmaan::magmaan(spec_h0, stats, estimator = magmaan_estimator)
    if (!isTRUE(mag_h1$converged) || !isTRUE(mag_h0$converged)) return(out)

    mag_single <- magmaan::fmg_tests_ordinal(
      mag_h1, stats, tests = single_tests, weight = magmaan_estimator
    )
    out["single"] <- p_diff(
      semTests::pvalues(lav_h1, single_tests), mag_single$p_value
    )
    mag_nested <- magmaan::fmg_nested_ordinal(
      mag_h1, mag_h0, stats, tests = nested_tests, weight = magmaan_estimator,
      A.method = "delta"
    )
    out["nested"] <- p_diff(
      semTests::pvalues_nested(
        lav_h0, lav_h1, tests = nested_tests, method = "2000",
        A.method = "delta"
      ),
      mag_nested$p_value
    )
  }, silent = TRUE)
  out
}

run_multigroup_categorical <- function(seed) {
  out <- na_pair
  try({
    data <- ordinalize(continuous_draw(seed))
    data$group <- rep(c("g1", "g2"), length.out = nrow(data))
    labels <- unique(as.character(data$group))
    lav_h1 <- lavaan::cfa(
      h1_model, data, ordered = ordered_names, group = "group",
      estimator = "WLSMV"
    )
    lav_h0 <- lavaan::cfa(
      h1_model, data, ordered = ordered_names, group = "group",
      group.equal = "loadings", estimator = "WLSMV"
    )
    if (!lavInspect(lav_h1, "converged") ||
        !lavInspect(lav_h0, "converged")) {
      return(out)
    }
    single_tests <- tests_for_df(lavaan::fitmeasures(lav_h1, "df"))
    nested_tests <- tests_for_df(nested_df_of(lav_h0, lav_h1))
    spec_h1 <- magmaan::model_spec(
      h1_model, ordered = ordered_names, group = "group",
      group_labels = labels
    )
    spec_h0 <- magmaan::model_spec(
      h1_model, ordered = ordered_names, group = "group",
      group_labels = labels, group_equal = "loadings"
    )
    stats <- magmaan::magmaan_core$data_ordinal_stats_from_df(
      data, spec_h1, ordered = ordered_names, group = "group",
      missing = "listwise", pd_gamma = "overlap", full_wls_weight = FALSE
    )
    mag_h1 <- magmaan::magmaan(spec_h1, stats, estimator = "DWLS")
    mag_h0 <- magmaan::magmaan(spec_h0, stats, estimator = "DWLS")
    if (!isTRUE(mag_h1$converged) || !isTRUE(mag_h0$converged)) return(out)

    mag_single <- magmaan::fmg_tests_ordinal(
      mag_h1, stats, tests = single_tests, weight = "DWLS"
    )
    out["single"] <- p_diff(
      semTests::pvalues(lav_h1, single_tests), mag_single$p_value
    )
    mag_nested <- magmaan::fmg_nested_ordinal(
      mag_h1, mag_h0, stats, tests = nested_tests, weight = "DWLS",
      A.method = "delta"
    )
    out["nested"] <- p_diff(
      semTests::pvalues_nested(
        lav_h0, lav_h1, tests = nested_tests, method = "2000",
        A.method = "delta"
      ),
      mag_nested$p_value
    )
  }, silent = TRUE)
  out
}

run_mixed_categorical <- function(seed) {
  out <- na_pair
  try({
    ord3 <- paste0("y", seq_len(3L))
    data <- ordinalize(continuous_draw(seed), vars = ord3)
    lav_h1 <- lavaan::cfa(
      h1_model, data, ordered = ord3, estimator = "WLSMV",
      meanstructure = TRUE
    )
    lav_h0 <- lavaan::cfa(
      h0_model, data, ordered = ord3, estimator = "WLSMV",
      meanstructure = TRUE
    )
    if (!lavInspect(lav_h1, "converged") ||
        !lavInspect(lav_h0, "converged")) {
      return(out)
    }
    single_tests <- tests_for_df(lavaan::fitmeasures(lav_h1, "df"))
    nested_tests <- tests_for_df(nested_df_of(lav_h0, lav_h1))
    spec_h1 <- magmaan::model_spec(
      h1_model, ordered = ord3, meanstructure = TRUE
    )
    spec_h0 <- magmaan::model_spec(
      h0_model, ordered = ord3, meanstructure = TRUE
    )
    stats <- magmaan::magmaan_core$data_mixed_ordinal_stats_from_df(
      data, spec_h1, ordered = ord3, missing = "listwise",
      full_wls_weight = FALSE
    )
    mag_h1 <- magmaan::magmaan(spec_h1, stats, estimator = "DWLS")
    mag_h0 <- magmaan::magmaan(spec_h0, stats, estimator = "DWLS")
    if (!isTRUE(mag_h1$converged) || !isTRUE(mag_h0$converged)) return(out)

    mag_single <- magmaan::fmg_tests_mixed_ordinal(
      mag_h1, stats, tests = single_tests, weight = "DWLS"
    )
    out["single"] <- p_diff(
      semTests::pvalues(lav_h1, single_tests), mag_single$p_value
    )
    mag_nested <- magmaan::fmg_nested_mixed_ordinal(
      mag_h1, mag_h0, stats, tests = nested_tests, weight = "DWLS",
      A.method = "delta"
    )
    out["nested"] <- p_diff(
      semTests::pvalues_nested(
        lav_h0, lav_h1, tests = nested_tests, method = "2000",
        A.method = "delta"
      ),
      mag_nested$p_value
    )
  }, silent = TRUE)
  out
}

# ---- cell registry: name, seed offset, tolerance, runner ---------------------

cells <- list(
  list(name = "continuous_GLS", offset = 1L, tol = 2e-3,
       run = function(s) run_continuous_ls(s, "GLS", "empirical")),
  list(name = "continuous_ULS", offset = 6L, tol = 2e-3,
       run = function(s) run_continuous_ls(s, "ULS", "normal")),
  list(name = "fiml_observed", offset = 2L, tol = 2e-3,
       run = function(s) run_fiml(s)),
  list(name = "categorical_DWLS_complete", offset = 3L, tol = 2e-3,
       run = function(s) run_categorical(s, pairwise = FALSE)),
  list(name = "categorical_DWLS_pairwise", offset = 4L, tol = 2e-3,
       run = function(s) run_categorical(s, pairwise = TRUE)),
  list(name = "categorical_ULS_complete", offset = 8L, tol = 2e-3,
       run = function(s) run_categorical(s, FALSE, "ULSMV", "ULS")),
  list(name = "mixed_categorical_DWLS", offset = 10L, tol = 2e-3,
       run = function(s) run_mixed_categorical(s)),
  list(name = "multigroup_categorical_DWLS", offset = 5L, tol = 5e-3,
       run = function(s) run_multigroup_categorical(s))
)

if (!identical(cells_filter, "all")) {
  keep <- strsplit(cells_filter, ",")[[1]]
  cells <- Filter(function(cell) cell$name %in% keep, cells)
}

cat(sprintf(
  "semTests %s | lavaan %s | magmaan %s\n",
  as.character(utils::packageVersion("semTests")),
  as.character(utils::packageVersion("lavaan")),
  as.character(utils::packageVersion("magmaan"))
))
cat(sprintf(
  "reps = %d | N = %d | cores = %d | tests = SB, SS, ALL, PEBA(min(4,df))\n\n",
  reps, n, ncores
))

seed_base <- 202607200L

summary_rows <- lapply(cells, function(cell) {
  draws <- parallel::mclapply(
    seq_len(reps),
    function(i) cell$run(seed_base + cell$offset * 100000L + i),
    mc.cores = ncores, mc.preschedule = TRUE
  )
  mat <- do.call(rbind, draws)
  worst <- apply(mat, 2L, max, na.rm = TRUE)
  worst[!is.finite(worst)] <- NA_real_
  pass <- vapply(c("single", "nested"), function(kind) {
    mean(mat[, kind] < cell$tol, na.rm = TRUE)
  }, numeric(1))
  row <- data.frame(
    cell = cell$name,
    tolerance = cell$tol,
    valid_single = sum(is.finite(mat[, "single"])),
    valid_nested = sum(is.finite(mat[, "nested"])),
    worst_single = worst[["single"]],
    worst_nested = worst[["nested"]],
    pass_rate_single = pass[["single"]],
    pass_rate_nested = pass[["nested"]],
    stringsAsFactors = FALSE
  )
  cat(sprintf(
    "%-30s worst single=%.2e nested=%.2e | pass s=%.0f%% n=%.0f%% | valid %d/%d\n",
    cell$name, row$worst_single, row$worst_nested,
    100 * row$pass_rate_single, 100 * row$pass_rate_nested,
    row$valid_single, reps
  ))
  row
})

summary_table <- do.call(rbind, summary_rows)
summary_table$reps <- reps
summary_table$sample_n <- n
summary_table$lavaan_version <- as.character(utils::packageVersion("lavaan"))
summary_table$magmaan_version <- as.character(utils::packageVersion("magmaan"))

if (nzchar(sim_output)) {
  dir.create(dirname(sim_output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary_table, sim_output, row.names = FALSE)
  cat("\nMachine-readable parity-simulation report:", sim_output, "\n")
}

worst_overall <- max(
  c(summary_table$worst_single, summary_table$worst_nested), na.rm = TRUE
)
cat(sprintf(
  "\nWorst parity difference across all cells: %.2e\n", worst_overall
))
