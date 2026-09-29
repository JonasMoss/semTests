# semTests (development version)

## New features

* `semtest()` builds score and criterion tests, for overall fit and for nested
  restrictions, as inspectable `semTests_test` objects. `pvalues()` now
  dispatches on these as well as on lavaan fits, preserving existing lavaan
  defaults and output names. Score tests accept an alternative fitted with
  `do.fit = FALSE`, since only the restricted fit has to be optimized.

* `semtest()`'s `weight` argument selects the matrix in the quadratic form, and
  therefore the statistic, separately from the curvature used to project out
  nuisance directions, which enters the reference law. `weight =
  "score_covariance"` studentizes fully and gives a chi-square reference; the
  curvature weightings give a weighted chi-square, so the usual calibrations
  apply. `weight = "expected_curvature"` with `curvature = "observed"` keeps the
  observed sensitivity the limit requires under misspecification while weighting
  by the more stable expected matrix.

* `inference_control()` collects the convention choices: `curvature`,
  `score_covariance`, `weight_effect`, and the numerical tolerances. Unset
  conventions follow the null being tested rather than one global default.
  Overall fit assumes the model holds, so it uses expected curvature with
  uncentered score products; nested restriction tests allow a misspecified
  baseline and use observed curvature with centered products. The resolved
  settings are recorded in the result's provenance.

* Score extraction covers ML and FIML likelihood scores including fitted-mean
  residuals, with missing-data contributions handled by pattern, and continuous
  or all-ordered least squares. Estimated first-stage weights and their joint
  influence are propagated by default; `weight_effect = "fixed"` treats the
  weights as known. Ordinal inference requires complete data and nonempty
  pairwise cells.

* `score_covariance = "expected"` puts model-based information in the covariance
  slot. Under complete data the information equality makes this the expected
  curvature, so the test reduces to `tests = "STD"`; it is a distinct test only
  under missing data, where the expectation also runs over the missingness
  patterns.

* A new vignette, `vignette("score-criterion", package = "semTests")`, explains
  the score and criterion tests, their hypotheses, and their settings.

## Deprecated

* `fiml.convention` is deprecated and will be removed in a future release.
  Supplying it, with either value, now gives a warning of class
  `semTests_warning_deprecated`; results are unchanged. Its `"lavaan"` option
  read lavaan's robust-test matrices as they were, and those change with the
  fit's `information` and `h1.information` options, so the same call on the
  same data could return different numbers depending on how the model was
  fitted. To reproduce lavaan's own tests, use `lavaan::lavTest()` or
  `lavaan::lavTestLRT()`.

## Performance

* Requests for several p-values share statistics and eigenvalue spectra within
  each call. Standard-only requests skip spectrum construction. Test ordering,
  duplicate requests, and calibration formulas are preserved.
* Expected-information single-model spectra use a symmetric residual-space
  factorization with groupwise covariance products. The observed-convention
  FIML spectrum uses the same reduction when its saturated metric is positive
  definite. Nested spectra use Cholesky solves when their metric permits it.
  Existing full-matrix constructions remain available as references and for
  conventions that do not admit these factorizations.
* Complete-data normal-likelihood ML can project empirical moment influences
  in row tiles when Gamma is not already cached, avoiding its full construction.
  Cached or user-supplied Gamma is reused. Biased nested covariance products
  are accumulated by group without a global block-diagonal Gamma matrix.

## Bug fixes

* Nested comparisons involving a nonlinear equality constraint used a constraint
  Jacobian evaluated at the starting values rather than the estimates, so the
  restriction was taken at the wrong point. Linear constraints were unaffected,
  their Jacobian being constant, and continuous data was protected by lavaan
  refusing nonlinear equality constraints upstream; categorical nested tests
  were affected.

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
