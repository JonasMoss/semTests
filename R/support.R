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
#' | DWLS family      | ordered/mixed | single or multi | listwise/pairwise | experimental |
#' | ULS family       | ordered/mixed | single or multi | listwise/pairwise | experimental |
#' | WLS (ADF)        | continuous  | single or multi | complete | degenerate (no-op)|
#' | WLS (ADF)        | ordered/mixed | single or multi | listwise/pairwise | degenerate (no-op)|
#'
#' ## Nested (`pvalues_nested()`)
#'
#' | Estimators              | Data        | Groups          | Missing          | Method      | A.method        | Status    |
#' |-------------------------|-------------|-----------------|------------------|-------------|-----------------|-----------|
#' | ML/MLM/MLR, GLS, ULS    | continuous  | single or multi | complete         | 2000 or 2001| --              | supported |
#' | ML / MLR (FIML)         | continuous  | single only     | FIML (both fits) | 2000 only   | exact or delta  | supported |
#' | DWLS/ULS/WLS families   | ordered/mixed | single or multi | listwise/pairwise | 2000 only | delta only | experimental |
#' | any                     | continuous  | --              | mixed / non-FIML | --          | --              | rejected  |
#'
#' ## Stability
#'
#' The classical normal-theory ML path (continuous, complete data) is the mature,
#' paper-backed core and is stable. Support for the other estimators (GLS, ULS,
#' categorical DWLS/ULS/WLS fits, FIML missing-data fits, and nested comparison
#' under FIML or categorical estimation) is **experimental** as of 0.9.0: the
#' methodology rests on the
#' references but the implementation surface is newer and less Monte-Carlo-vetted,
#' so the API and numerical details may change. Everything in the tables above is
#' validated and tested; configurations outside them are refused at the entry
#' point.
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
#' * **Categorical and mixed-indicator fits** use lavaan's biased inspected
#'   `UGamma` spectrum and unscaled estimator statistic directly. Supported
#'   estimator families are DWLS (including WLSMV/WLSM/WLSMVS), ULS (including
#'   ULSMV), and full WLS. Nested comparison is restricted to
#'   `method = "2000"` and `A.method = "delta"`.
#' * **Pairwise categorical missingness** means support for lavaan's pairwise
#'   sample-statistic and `UGamma` calculation. It is not FIML and is not a
#'   claim of generally MAR-valid inference.
#' * **FIML** (missing data) uses the biased gamma and the standard statistic
#'   only; `UG` and `RLS` are refused. FIML requires a single group, continuous
#'   data, and no fixed exogenous covariates. `fiml.convention = "observed"`
#'   (the default) uses observed saturated and model information, independently
#'   validated against magmaan. `"lavaan"` reproduces lavaan 0.7-2's inspected
#'   robust-test spectrum, including its expected H1 weight.
#'
#' ## Why ADF/WLS is degenerate
#'
#' For full WLS (ADF), with continuous or categorical sample statistics, the
#' model test statistic already uses the
#' asymptotically-correct weight matrix, so its null distribution is the exact
#' chi-square and the eigenvalue correction collapses to the identity: every
#' p-value equals the ordinary `1 - pchisq(chisq, df)`. It is accepted but adds
#' nothing.
#'
#' @seealso [pvalues()], [pvalues_nested()]
#' @name semTests-support
#' @keywords internal
NULL
