#!/usr/bin/env Rscript

# Plan item B4: structured vs unstructured saturated FIML information in the
# goodness-of-fit eigenvalue tests. Seeded, non-CRAN, resumable.
#
# The released pvalues() path for FIML fits (fiml_lambdas(), "observed"
# convention) evaluates the observed saturated information at the unrestricted
# EM moments (lavaan h1.information = "unstructured"). semtest() evaluates it at
# the model-implied moments, which is lavaan's h1.information = "structured"
# (identical, entry for entry, to gaussian_pattern_parts()$observed / n).
# Only that matrix differs between the two conventions here: Gamma, Delta and
# the statistic are shared, and both are computed on the same draw, so the
# comparison is paired.
#
# Design (12 data cells x 2 conventions = 24 cells):
#   N in {100, 200, 400}; normal or t5 (radial, variance-preserving);
#   correct model, or a population with an omitted .20 cross-loading
#   (y4 on f1); MAR missingness as in tools/fiml-validation-simulation.R.
#   Tests SB, pEBA4, PALL; rejection at .01, .05, .10.
#
# Run from the package root:
#   REPS=2000 CHUNK=250 NCORES=4 OUT_DIR=workspace/fiml-h1-information \
#     nice -n 10 Rscript tools/fiml-h1-information-simulation.R
# Each finished chunk is one file under OUT_DIR/chunks, written atomically, so
# an interrupted run resumes where it stopped. The summary is written once all
# chunks exist.

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to run against the current semTests source tree.",
       call. = FALSE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(lavaan)
  library(parallel)
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
})

reps <- as.integer(Sys.getenv("REPS", "2000"))
chunk_size <- as.integer(Sys.getenv("CHUNK", "250"))
ncores <- as.integer(Sys.getenv("NCORES", "4"))
out_dir <- Sys.getenv("OUT_DIR", "workspace/fiml-h1-information")
chunk_dir <- file.path(out_dir, "chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

population_correct <- "
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
population_misspecified <- sub(
  "f2 =~ .7*y4", "f1 =~ .2*y4\n  f2 =~ .7*y4", population_correct,
  fixed = TRUE
)
fit_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
"

cells <- expand.grid(
  n = c(100L, 200L, 400L),
  distribution = c("normal", "t5"),
  model = c("correct", "misspecified"),
  stringsAsFactors = FALSE
)
cells$cell <- seq_len(nrow(cells))

conventions <- c("unstructured", "structured")
tests <- c("sb", "peba4", "pall")
value_names <- c(
  "converged", "chisq",
  as.vector(outer(conventions, tests, paste, sep = "_"))
)

make_data <- function(cell, rep) {
  set.seed(2026102800L + 100000L * cell$cell + rep)
  population <- if (cell$model == "correct") {
    population_correct
  } else {
    population_misspecified
  }
  data <- simulateData(population, sample.nobs = cell$n, model.type = "cfa",
                       std.lv = TRUE)
  if (cell$distribution == "t5") {
    radial <- sqrt(3 / stats::rchisq(cell$n, df = 5))
    data[] <- lapply(data, function(x) x * radial)
  }
  z <- as.numeric(scale(data$y6))
  miss1 <- stats::runif(cell$n) < stats::plogis(stats::qlogis(.30) + .65 * z)
  miss2 <- stats::runif(cell$n) < stats::plogis(stats::qlogis(.30) - .65 * z)
  data[miss1, c("y1", "y2", "y3")] <- NA_real_
  data[miss2, c("y4", "y5")] <- NA_real_
  data
}

# Observed saturated information under either lavaan h1.information setting,
# on the scale fiml_h1_information_observed() uses.
h1_information <- function(fit, convention) {
  object <- fit
  object@Options$h1.information[] <- convention
  semTests:::fiml_group_bdiag(
    lavaan::lavTech(object, "h1.information.observed"),
    fit, "observed H1 information", scale = "weight"
  )
}

# The "observed" branch of semTests:::fiml_lambdas() with H supplied.
lambdas_with <- function(fit, df, H) {
  Gamma <- semTests:::fiml_gamma_matrix(fit)
  Delta <- semTests:::fiml_delta_effective(fit)
  B <- semTests:::residual_metric_factor(list(H), Delta, df)
  if (!is.null(B)) {
    return(semTests:::positive_symmetric_eigenvalues(
      crossprod(B, Gamma %*% B), df, "FIML residual spectrum"
    ))
  }
  HD <- H %*% Delta
  P <- semTests:::symmetrize(crossprod(Delta, HD))
  U <- H - HD %*% semTests:::generalized_inverse(P) %*% t(HD)
  semTests:::fiml_sandwich_eigenvalues(U, Gamma, df)
}

one <- function(cell, rep) {
  out <- stats::setNames(rep(NA_real_, length(value_names)), value_names)
  out["converged"] <- 0
  data <- make_data(cell, rep)
  fit <- tryCatch(
    suppressWarnings(cfa(fit_model, data, missing = "fiml",
                         estimator = "MLR", std.lv = TRUE)),
    error = function(e) NULL
  )
  if (is.null(fit) || !lavInspect(fit, "converged")) return(out)
  out["converged"] <- 1
  df <- as.integer(fitMeasures(fit, "df"))
  chisq <- as.numeric(fitMeasures(fit, "chisq"))
  out["chisq"] <- chisq
  for (convention in conventions) {
    lambda <- tryCatch(
      lambdas_with(fit, df, h1_information(fit, convention)),
      error = function(e) NULL
    )
    if (is.null(lambda)) next
    out[paste0(convention, "_sb")] <-
      semTests:::trad_pvalue(df, chisq, lambda, "sb")
    out[paste0(convention, "_pall")] <-
      semTests:::trad_pvalue(df, chisq, lambda, "pall")
    out[paste0(convention, "_peba4")] <-
      semTests:::peba_pvalue(chisq, lambda, 4L)
  }
  out
}

# Self-check: the unstructured branch must reproduce the released pvalues().
local({
  cell <- cells[cells$n == 400L & cells$distribution == "normal" &
                  cells$model == "correct", ]
  data <- make_data(cell, 0L)
  fit <- cfa(fit_model, data, missing = "fiml", estimator = "MLR",
             std.lv = TRUE)
  mine <- one(cell, 0L)[paste0("unstructured_", tests)]
  released <- suppressWarnings(as.numeric(
    semTests::pvalues(fit, c("SB", "pEBA4", "PALL"))
  ))
  gap <- max(abs(mine - released))
  if (!is.finite(gap) || gap > 1e-10) {
    stop(sprintf("Unstructured branch differs from pvalues() by %.3g.", gap),
         call. = FALSE)
  }
  cat(sprintf("Self-check against pvalues(): max gap %.2e\n", gap))
})

chunks <- expand.grid(
  chunk = seq_len(ceiling(reps / chunk_size)), cell = cells$cell
)
chunk_file <- function(cell, chunk) {
  file.path(chunk_dir, sprintf("cell-%02d-chunk-%02d.csv", cell, chunk))
}
todo <- which(!file.exists(mapply(chunk_file, chunks$cell, chunks$chunk)))
cat(sprintf(
  "semTests %s | lavaan %s | %d reps/cell in chunks of %d | %d of %d chunks to run on %d cores\n",
  as.character(utils::packageVersion("semTests")),
  as.character(utils::packageVersion("lavaan")),
  reps, chunk_size, length(todo), nrow(chunks), ncores
))

run_chunk <- function(k) {
  cell <- cells[chunks$cell[k], ]
  first <- (chunks$chunk[k] - 1L) * chunk_size + 1L
  ids <- seq.int(first, min(reps, first + chunk_size - 1L))
  values <- do.call(rbind, lapply(ids, function(r) one(cell, r)))
  frame <- data.frame(cell = cell$cell, rep = ids, values, row.names = NULL)
  target <- chunk_file(cell$cell, chunks$chunk[k])
  temporary <- paste0(target, ".tmp")
  utils::write.csv(frame, temporary, row.names = FALSE)
  file.rename(temporary, target)
  invisible(target)
}

started <- Sys.time()
status <- parallel::mclapply(
  todo, function(k) tryCatch(run_chunk(k), error = function(e) conditionMessage(e)),
  mc.cores = ncores, mc.preschedule = FALSE
)
failed <- !vapply(status, function(s) is.character(s) && file.exists(s), logical(1))
if (any(failed)) {
  cat("Failed chunks:\n")
  print(unlist(status[failed]))
}
cat(sprintf("Ran %d chunks in %.1f minutes\n", length(todo),
            as.numeric(difftime(Sys.time(), started, units = "mins"))))

files <- mapply(chunk_file, chunks$cell, chunks$chunk)
if (!all(file.exists(files))) {
  cat("Not all chunks are present yet; rerun to resume.\n")
  quit(status = 1L)
}

raw <- do.call(rbind, lapply(files, utils::read.csv))
raw <- merge(cells, raw, by = "cell")
utils::write.csv(raw, file.path(out_dir, "raw.csv"), row.names = FALSE)

summary_rows <- list()
for (cell_id in cells$cell) {
  d <- raw[raw$cell == cell_id, ]
  for (test in tests) {
    pu <- d[[paste0("unstructured_", test)]]
    ps <- d[[paste0("structured_", test)]]
    both <- is.finite(pu) & is.finite(ps)
    for (alpha in c(.01, .05, .10)) {
      ru <- pu[both] < alpha
      rs <- ps[both] < alpha
      diff <- rs - ru
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        cells[cells$cell == cell_id, c("n", "distribution", "model")],
        test = test, alpha = alpha,
        reps = nrow(d), converged = sum(d$converged == 1),
        valid_pairs = sum(both),
        reject_unstructured = mean(ru),
        reject_structured = mean(rs),
        mcse_unstructured = sqrt(mean(ru) * (1 - mean(ru)) / sum(both)),
        mcse_structured = sqrt(mean(rs) * (1 - mean(rs)) / sum(both)),
        paired_diff = mean(diff),
        paired_mcse = stats::sd(diff) / sqrt(sum(both)),
        discordant = sum(diff != 0),
        max_abs_p_gap = max(abs(ps[both] - pu[both])),
        row.names = NULL
      )
    }
  }
}
summary_table <- do.call(rbind, summary_rows)
summary_table$lavaan_version <- as.character(utils::packageVersion("lavaan"))
utils::write.csv(summary_table, file.path(out_dir, "summary.csv"),
                 row.names = FALSE)
cat("Summary written to", file.path(out_dir, "summary.csv"), "\n")
print(summary_table[summary_table$alpha == .05,
                    c("n", "distribution", "model", "test", "valid_pairs",
                      "reject_unstructured", "reject_structured",
                      "paired_diff", "paired_mcse")],
      row.names = FALSE, digits = 3)
