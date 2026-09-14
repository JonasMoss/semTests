## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(semTests)

## ----models, echo=FALSE, out.width="100%", fig.alt="Two one-factor models with four indicators. The unrestricted model has loadings a and b on x2 and x3; the restricted model uses the same loading on both paths.", fig.cap="The same factor model before and after imposing equality. Blue paths identify the restriction. Residuals are uncorrelated; intercepts and variances are omitted from the drawing."----
knitr::include_graphics(if (knitr::is_latex_output()) {
  "figures/score-wald-models.pdf"
} else {
  "figures/score-wald-models.svg"
})

## ----fits---------------------------------------------------------------------
d <- lavaan::HolzingerSwineford1939
model <- "f =~ x1 + a*x2 + b*x3 + x4"
fit <- lavaan::cfa(model, d, meanstructure = TRUE)
restricted <- lavaan::cfa(paste(model, "a == b", sep = "\n"),
                          d, meanstructure = TRUE)
lavaan::coef(fit)[c("a", "b")]

## ----wald---------------------------------------------------------------------
R <- matrix(0, 1, length(lavaan::coef(fit)))
R[1, match(c("a", "b"), names(lavaan::coef(fit)))] <- c(1, -1)
w <- wald_test(fit, R, vcov = "sandwich",
               curvature = "observed", score_covariance = "centered",
               hypothesis = "a = b", keep = TRUE)
w

## ----score--------------------------------------------------------------------
profile <- score_test_nested(restricted, fit,
                             metric = "sandwich",
                             curvature = "observed", score_covariance = "centered",
                             hypothesis = "a = b", keep = TRUE)
profile

## ----template, eval=FALSE-----------------------------------------------------
# template <- lavaan::cfa(model, d, meanstructure = TRUE, do.fit = FALSE)
# score_test_nested(restricted, template, metric = "sandwich",
#                   curvature = "observed", score_covariance = "centered")

## ----criterion----------------------------------------------------------------
criterion <- criterion_test_nested(restricted, fit,
                                    curvature = "observed", score_covariance = "centered",
                                    hypothesis = "a = b")
criterion

## ----comparison, echo=FALSE---------------------------------------------------
comparison <- data.frame(
  Test = c("Wald", "Score", "Criterion difference"),
  Statistic = c(w$statistic, profile$statistic, criterion$statistic),
  `p-value` = c(w$p_value, profile$p_value, criterion$p_value),
  check.names = FALSE
)
knitr::kable(comparison, digits = 4)

## ----calibrations-------------------------------------------------------------
curvature_score <- score_test_nested(restricted, fit,
                                     metric = "information",
                                     curvature = "observed", score_covariance = "centered")
pvalues(curvature_score, c("ALL", "SB"))

## ----global-------------------------------------------------------------------
global <- score_test(fit, metric = "information",
                      curvature = "expected", score_covariance = "uncentered")
global

## ----weights------------------------------------------------------------------
fit_dwls <- lavaan::cfa(model, d, meanstructure = TRUE, estimator = "DWLS",
                         ordered = FALSE)
restricted_dwls <- lavaan::cfa(paste(model, "a == b", sep = "\n"), d,
                              meanstructure = TRUE, estimator = "DWLS",
                         ordered = FALSE)
weighted_score <- score_test_nested(restricted_dwls, fit_dwls,
                                    metric = "sandwich",
                                    curvature = "observed", score_covariance = "centered",
                                    weight_effect = "estimated")
weighted_score

