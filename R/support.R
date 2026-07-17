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
#' | Estimator        | Data        | Groups          | Missing  | Availability | Stability |
#' |------------------|-------------|-----------------|----------|--------------|-----------|
#' | ML / MLM / MLR   | continuous  | single or multi | complete | available    | stable |
#' | GLS              | continuous  | single or multi | complete | available    | experimental |
#' | ULS              | continuous  | single or multi | complete | available    | experimental |
#' | ML / MLR (FIML)  | continuous  | single only     | FIML     | available    | experimental |
#' | DWLS family      | ordered/mixed | single or multi | listwise/pairwise | available | experimental |
#' | ULS family       | ordered/mixed | single or multi | listwise/pairwise | available | experimental |
#' | WLS (ADF)        | continuous  | single or multi | complete | available (no-op) | experimental |
#' | WLS (ADF)        | ordered/mixed | single or multi | listwise/pairwise | available (no-op) | experimental |
#'
#' ## Nested (`pvalues_nested()`)
#'
#' | Estimators              | Data        | Groups          | Missing          | Method      | A.method        | Availability | Stability |
#' |-------------------------|-------------|-----------------|------------------|-------------|-----------------|--------------|-----------|
#' | ML/MLM/MLR              | continuous  | single or multi | complete         | 2000 or 2001| --              | available | stable |
#' | GLS, ULS                | continuous  | single or multi | complete         | 2000 or 2001| --              | available | experimental |
#' | ML / MLR (FIML)         | continuous  | single only     | FIML (both fits) | 2000 only   | exact or delta  | available | experimental |
#' | DWLS/ULS/WLS families   | ordered/mixed | single or multi | listwise/pairwise | 2000 only | delta only | available | experimental |
#' | any                     | continuous  | --              | mixed / non-FIML | --          | --              | rejected | -- |
#'
#' ## Stability
#'
#' The classical normal-theory ML path (continuous, complete data) is the mature,
#' paper-backed core and is stable. Support for the other estimators (GLS, ULS,
#' categorical DWLS/ULS/WLS fits, FIML missing-data fits, and nested comparison
#' under FIML or categorical estimation) is **experimental** as of 0.9.0: the
#' methodology rests on the
#' references but the implementation surface is newer and less Monte-Carlo-vetted,
#' so the API and numerical details may change. Available rows have deterministic
#' implementation tests; the experimental rows have not received the same breadth
#' of Monte Carlo evaluation as the paper-backed classical ML core. Configurations
#' outside the tables are refused at the entry point.
#'
#' ## Fit quality and nested comparability
#'
#' Both entry points require converged lavaan fits and warn when lavaan's
#' post-estimation check reports an inadmissible solution. Nested comparisons
#' additionally require the same requested estimator, fitting conventions,
#' observed variables, groups, sample sizes, raw data, and missingness mask.
#' These checks establish comparability but cannot prove that two substantively
#' different model specifications are genuinely nested.
#'
#' The constrained model belongs in `m0`. If the models are supplied in reverse
#' df order, semTests warns and swaps them. For `method = "2001"`, materially
#' negative leading eigenvalues produce an error with a recommendation to inspect
#' the fits and use method 2000; semTests never silently changes the method.
#'
#' ## Statistic and gamma options
#'
#' Test names use one of `TEST`, `TEST_UG`, `TEST_ML`, `TEST_RLS`,
#' `TEST_UG_ML`, or `TEST_UG_RLS` (e.g. `"SB_UG_RLS"`); see [pvalues()] for the
#' test families. Two of the options are defined only for the classical case:
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
#'   (the default) uses observed saturated and model information, following the
#'   recommendation of Kenward and Molenberghs (1998), and is independently
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
#' @references
#'
#' Foldnes, N., Moss, J., & Grønneberg, S. (2025). Improved goodness of fit
#' procedures for structural equation models. *Structural Equation Modeling: A
#' Multidisciplinary Journal*, 32(1), 1--13.
#' \doi{10.1080/10705511.2024.2372028}
#'
#' Foldnes, N., Grønneberg, S., & Moss, J. (2026). Penalized eigenvalue block
#' averaging: Extension to nested model comparison and Monte Carlo evaluations.
#' *Behavior Research Methods*, 58, article 107.
#' \doi{10.3758/s13428-026-02968-4}
#'
#' Kenward, M. G., & Molenberghs, G. (1998). Likelihood based frequentist
#' inference when data are missing at random. *Statistical Science*, 13(3),
#' 236--247. \doi{10.1214/ss/1028905886}
#'
#' Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
#' analysis of moment structures. In R. D. H. Heijmans, D. S. G. Pollock, &
#' A. Satorra (Eds.), *Innovations in Multivariate Statistical Analysis*
#' (pp. 233--247). Kluwer Academic.
#' \doi{10.1007/978-1-4615-4603-0_17}
#'
#' @seealso [pvalues()], [pvalues_nested()]
#' @name semTests-support
#' @keywords models
NULL
