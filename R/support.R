#' Supported estimators, data types, and configurations
#'
#' Use this page to check what [pvalues()] and [pvalues_nested()] support.
#' The eigenvalue-based p-values target a limiting null law that is a weighted
#' sum of chi-squares. This result holds for any minimum-discrepancy estimator
#' and therefore includes estimators beyond normal-theory ML. The tables below
#' give the supported combinations. Configurations outside them are refused at
#' the entry point (see [check_supported()]), and the test suite follows the
#' same boundary.
#'
#' @details
#'
#' ## Single-model (`pvalues()`)
#'
#' | Estimator        | Data        | Groups          | Missing  | Availability | Stability |
#' |------------------|-------------|-----------------|----------|--------------|-----------|
#' | ML / MLM / MLR   | continuous  | single or multi | complete | available    | stable |
#' | GLS              | continuous  | single or multi | complete | available    | stable |
#' | ULS              | continuous  | single or multi | complete | available    | stable |
#' | ML / MLR (FIML)  | continuous  | single or multi | FIML     | available    | stable |
#' | DWLS family      | ordered/mixed | single or multi | listwise/pairwise | available | stable |
#' | ULS family       | ordered/mixed | single or multi | listwise/pairwise | available | stable |
#' | WLS (ADF)        | continuous  | single or multi | complete | rejected | -- |
#' | WLS (ADF)        | ordered/mixed | single or multi | listwise/pairwise | rejected | -- |
#'
#' ## Nested (`pvalues_nested()`)
#'
#' | Estimators              | Data        | Groups          | Missing          | Method      | A.method        | Availability | Stability |
#' |-------------------------|-------------|-----------------|------------------|-------------|-----------------|--------------|-----------|
#' | ML/MLM/MLR              | continuous  | single or multi | complete         | 2000        | --              | available | stable |
#' | GLS, ULS                | continuous  | single or multi | complete         | 2000        | --              | available | stable |
#' | ML / MLR (FIML)         | continuous  | single or multi | FIML (both fits) | 2000 only   | exact or delta  | available | stable |
#' | DWLS/ULS families       | ordered/mixed | single or multi | listwise/pairwise | 2000 only | delta only | available | stable |
#' | any                     | continuous  | --              | mixed / non-FIML | --          | --              | rejected | -- |
#'
#' ## Stability
#'
#' The classical normal-theory ML path (continuous, complete data) is the mature,
#' paper-backed core. The other supported estimators -- GLS, ULS, categorical
#' DWLS/ULS fits, FIML missing-data fits, and their nested comparisons -- are
#' also marked stable: each is validated to numerical tolerance against an
#' independent implementation (magmaan) across a simulation battery, in addition
#' to the deterministic implementation tests every available row carries. The
#' `fiml.convention = "lavaan"` option reproduces lavaan's own spectrum for
#' compatibility rather than as an independently validated construction.
#' Configurations outside the tables are refused at the entry point.
#'
#' ## Observed exogenous covariates
#'
#' Observed exogenous predictors are supported when they are modeled jointly
#' with the other variables. In lavaan, request this random-x analysis with
#' `fixed.x = FALSE` and `conditional.x = FALSE`. The independent validation
#' suite covers continuous ML, GLS, ULS, FIML, and mixed categorical DWLS in
#' this setting, for single models and nested comparisons.
#'
#' Fixed or conditional observed exogenous predictors are currently refused.
#' Their reference distribution conditions on the realized covariate design
#' and needs a separate saturated-model projection. If a covariate is intended
#' to be fixed, keep that scientific choice and do not silently change it to
#' random merely to pass the software check.
#'
#' ## Fit quality and nested comparability
#'
#' Both entry points require converged lavaan fits and warn when lavaan's
#' post-estimation check reports an inadmissible solution. Nested comparisons
#' additionally require the same requested estimator, fitting conventions,
#' observed variables, groups, sample sizes, raw data, and missingness mask.
#' These checks establish comparability. The scientific argument that one model
#' is genuinely nested within the other still belongs in the analysis.
#'
#' The constrained model belongs in `m0`. If the models are supplied in reverse
#' df order, semTests warns and swaps them. Only the Satorra-2000 reduction is
#' available. `method = "2001"` has been withdrawn because it performs poorly.
#'
#' ## Statistic and gamma options
#'
#' Test names use one of `TEST`, `TEST_UG`, `TEST_ML`, `TEST_RLS`,
#' `TEST_UG_ML`, or `TEST_UG_RLS` (e.g. `"SB_UG_RLS"`). See [pvalues()] for the
#' test families. Two of the options are defined only for the classical case:
#'
#' * The **RLS** statistic (`browne.residual.nt.model`, Browne 1974) and the
#'   **`UG`** (Du-Bentler) unbiased gamma are available **only** for classical
#'   normal-theory ML with continuous, complete data, `estimator = "ML"`,
#'   `fixed.x = FALSE`, and `conditional.x = FALSE`. lavaan degrades RLS to ADF
#'   elsewhere, and the Du-Bentler correction has no derivation for the other
#'   families or a fixed covariate design. `semTests` refuses either request and
#'   points you to the standard statistic and biased gamma.
#' * **Categorical and mixed-indicator fits** use lavaan's biased inspected
#'   `UGamma` spectrum and unscaled estimator statistic directly. Supported
#'   estimator families are DWLS (including WLSMV/WLSM/WLSMVS) and ULS (including
#'   ULSMV). Nested comparison is restricted to `method = "2000"` and
#'   `A.method = "delta"`.
#' * **Pairwise categorical missingness** means support for lavaan's pairwise
#'   sample-statistic and `UGamma` calculation. Pairwise deletion has its own
#'   missingness assumptions and does not provide FIML or generally MAR-valid
#'   inference.
#' * **FIML** (missing data) uses the biased gamma and the standard statistic
#'   only. `UG` and `RLS` are refused. FIML supports one or several groups with
#'   continuous data. Observed exogenous predictors are supported under joint
#'   random-x inference, using `fixed.x = FALSE` and `conditional.x = FALSE`.
#'   `fiml.convention = "observed"` (the default) uses observed saturated and
#'   model information, following the recommendation of Kenward and Molenberghs
#'   (1998), and is independently validated against magmaan. `"lavaan"`
#'   reproduces lavaan 0.7-2's inspected robust-test spectrum, including its
#'   expected H1 weight.
#'
#' ## Why full WLS/ADF is refused
#'
#' For full WLS (ADF), with continuous or categorical sample statistics, the
#' fitting weight is already the inverse of the asymptotic moment covariance, so
#' `UGamma` is a projector: its nonzero eigenvalues are all exactly one and the
#' eigenvalue correction collapses to the identity. Every robust p-value would
#' therefore equal the ordinary `1 - pchisq(chisq, df)`. Because the correction
#' adds nothing, full WLS is refused at the entry point. Fit a DWLS/ULS family
#' for a non-degenerate robust test, or read the standard chi-square directly.
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
