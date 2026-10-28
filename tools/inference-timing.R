#!/usr/bin/env Rscript
# Optional timings against the pre-sharing implementation, with numerical parity.
# Run from the package root; REPS controls repetitions (default 10).
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
baseline <- "e54d90b27c03f151a73ad57a456308ad44862328"
old_source <- system2("git", c("show", paste0(baseline, ":R/pvalues.R")),
                      stdout = TRUE)
stopifnot(is.null(attr(old_source, "status")))
old <- new.env(parent = asNamespace("semTests"))
eval(parse(text = old_source), envir = old)
n <- as.integer(Sys.getenv("REPS", "10"))
stopifnot(length(n) == 1L, !is.na(n), n > 0L)
model <- "visual =~ x1+x2+x3
          textual =~ x4+x5+x6
          speed =~ x7+x8+x9"
fit <- lavaan::cfa(model, lavaan::HolzingerSwineford1939, estimator = "MLM")
elapsed <- function(f) {
  f()
  unname(system.time(for (i in seq_len(n)) f())[["elapsed"]]) * 1000 / n
}
rows <- lapply(list(
  standard = "STD_ML",
  shared_spectrum = c("SB_ML", "SS_ML", "SF_ML"),
  shared_rls = c("SB_RLS", "SS_RLS", "SF_RLS"),
  mixture = c("SB_ML", "ALL_ML", "PEBA4_ML")
), function(specs) {
  before <- function() old$pvalues(fit, specs)
  after <- function() semTests::pvalues(fit, specs)
  stopifnot(isTRUE(all.equal(before(), after(), tolerance = 1e-12)))
  data.frame(tests = paste(specs, collapse = ","),
             before_ms = elapsed(before), after_ms = elapsed(after))
})
result <- do.call(rbind, rows)
result$speedup <- result$before_ms / result$after_ms
dir.create("workspace", showWarnings = FALSE)
write.csv(result, "workspace/inference-timing.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "workspace/inference-timing-session.txt")
print(result)
cat("Baseline:", baseline, "; repetitions:", n, "\n")
