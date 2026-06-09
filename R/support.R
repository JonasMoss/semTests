#' Supported estimators, data types, and configurations
#'
#' The authoritative list of what [pvalues()] and [pvalues_nested()] support.
#' The eigenvalue-based p-values target the limiting null law of the test
#' statistic -- a weighted sum of chi-squares -- which holds for any
#' minimum-discrepancy estimator, so the tests are not specific to normal-theory
#' ML. This page is the single source of truth: anything not listed here is
#' refused at the entry point (see [check_supported()]), and the test suite is
#' written against it.
#'
#' @details
#'
#' ## Single-model (`pvalues()`)
#'
#' | Estimator        | Data        | Groups          | Missing  | Status            |
#' |------------------|-------------|-----------------|----------|-------------------|
#' | ML / MLM / MLR   | continuous  | single or multi | complete | supported         |
#' | GLS              | continuous  | single or multi | complete | supported         |
#' | ULS              | continuous  | single or multi | complete | supported         |
#' | ML / MLR (FIML)  | continuous  | single only     | FIML     | supported         |
#' | WLSMV / DWLS     | categorical | single or multi | complete | supported         |
#' | WLS (ADF)        | continuous  | single          | complete | degenerate (no-op)|
#'
#' ## Nested (`pvalues_nested()`)
#'
#' | Estimators              | Data        | Groups          | Missing          | Method      | A.method        | Status    |
#' |-------------------------|-------------|-----------------|------------------|-------------|-----------------|-----------|
#' | ML/MLM/MLR, GLS, ULS    | continuous  | single or multi | complete         | 2000 or 2001| --              | supported |
#' | ML / MLR (FIML)         | continuous  | single only     | FIML (both fits) | 2000 only   | exact or delta  | supported |
#' | WLSMV / DWLS            | categorical | any             | any              | --          | --              | rejected  |
#' | any                     | continuous  | --              | mixed / non-FIML | --          | --              | rejected  |
#'
#' ## Statistic and gamma options
#'
#' Test names have the form `(test)_(ug?)_(ml|rls?)` (e.g. `"SB_UG_RLS"`); see
#' [pvalues()] for the test families. Two of the options are defined only for the
#' classical case:
#'
#' * The **RLS** statistic (`browne.residual.nt.model`, Browne 1974) and the
#'   **`UG`** (Du-Bentler) unbiased gamma are available **only** for classical
#'   normal-theory ML: continuous, complete data, `estimator = "ML"`. Off that
#'   case they are undefined -- lavaan silently degrades RLS to ADF, and the
#'   Du-Bentler correction has no derivation -- so requesting them is refused.
#'   The standard statistic and the biased gamma are used instead.
#' * **FIML** (missing data) uses the biased gamma and the standard statistic
#'   only; `UG` and `RLS` are refused. The missing-data spectrum is renormalised
#'   (see [rescale_missing()]). FIML requires a single group, continuous data,
#'   and no fixed exogenous covariates, and -- because the inference rests on
#'   the observed information under MAR -- a fit with `information = "expected"`
#'   triggers a warning (see [warn_fiml_information()]).
#'
#' ## Why ADF/WLS is degenerate
#'
#' For full WLS (ADF) the model test statistic already uses the
#' asymptotically-correct weight matrix, so its null distribution is the exact
#' chi-square and the eigenvalue correction collapses to the identity: every
#' p-value equals the ordinary `1 - pchisq(chisq, df)`. It is accepted but adds
#' nothing.
#'
#' @seealso [pvalues()], [pvalues_nested()]
#' @name semTests-support
#' @keywords internal
NULL
