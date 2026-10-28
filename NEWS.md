# semTests 2.0.0

## Breaking changes

* `pvalues_nested()` has a new default restriction matrix,
  `A.method = "weighted"`, which does not depend on the units of the observed
  variables or on the factor identification. lavaan's construction remains
  available as `A.method = "delta"`, which reproduces
  `lavTestLRT(method = "satorra.2000")` but changes with the units.
* semTests requires `lavaan (>= 0.7-3)`. Its asymptotic covariance for mixed
  continuous and ordinal data is now robust for the continuous variables.

## New features

* `semtest()` computes score and criterion-difference tests of overall fit and
  of restrictions in nested models. The restriction tests remain valid when the
  larger model is misspecified. Overall fit covers complete-data ML and FIML;
  restriction tests also cover continuous least squares and complete,
  all-ordered DWLS, WLS, and ULS fits.
* `pvalues()` accepts `semtest()` results, so the usual calibrations (`ALL`,
  `SB`, `SS`, `PEBA<j>`, ...) apply to them.
* `inference_control()` collects the test conventions (curvature, score
  covariance, fitting-weight effect) and numerical settings.
* New vignette: `vignette("score-criterion", package = "semTests")`.

## Deprecated

* `fiml.convention` is deprecated. Its `"lavaan"` option depended on the fit's
  `information` and `h1.information` options; use `lavaan::lavTest()` or
  `lavaan::lavTestLRT()` to reproduce lavaan's tests.

## Bug fixes

* `pvalues_nested(A.method = "exact")` was ignored for complete-data fits.
* Nested comparisons with nonlinear equality constraints evaluated the
  constraint Jacobian at the starting values instead of the estimates. This
  affected categorical data.
* Nested FIML tests lost accuracy when the variables had very different
  scales.

## Other changes

* Requests for several tests share their statistics and eigenvalues, and
  complete-data ML no longer builds the full Gamma matrix.
* The vignettes are reorganized into four: an introduction, ordinal
  indicators, missing data, and score and criterion tests.

# semTests 1.0.0

Earlier versions marked several estimators as experimental. This release
promotes them to stable.

## New features

* GLS, ULS, categorical DWLS and ULS, mixed continuous and ordinal indicators,
  and FIML, each with one or several groups, now work for both single models and
  nested comparisons.
* A `fiml.convention` argument selects the FIML information convention. The
  default disagrees with lavaan but has better statistical properties; set
  `fiml.convention = "lavaan"` to reproduce lavaan's numbers exactly.
* New vignettes cover measurement invariance with ordinal indicators and latent
  growth models under FIML.

## Bug fixes

* The reweighted least-squares statistic was read as `NULL` under lavaan 0.7-2
  and is now correct. semTests requires `lavaan (>= 0.7-2)`.

## Removed

* Nested method 2001 is withdrawn because of its poor performance.
* Full WLS is refused because its correction reduces to the ordinary chi-square.
