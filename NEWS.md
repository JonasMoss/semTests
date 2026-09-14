# semTests 1.0.0

## Development performance improvements

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

Earlier versions marked several estimators as experimental. This release
promotes them to stable.

## New features

* New inference interfaces use `bread` for curvature and `meat` for score
  covariance, including their matrix inputs and returned ingredient/provenance
  fields. Defaults and numerical calculations are unchanged. Wald parameter
  covariance remains `vcov`, and score quadratic weighting remains `metric`.

* `score_test()` extracts ML/FIML likelihood scores, including fitted-mean
  residuals, and projects nuisance directions using expected or observed
  saturated curvature. Missing-data contributions are handled by pattern.
* `score_test_nested()` adds joint affine restriction releases at the restricted
  fit, with full criterion curvature. `criterion_test_nested()` adds the actual
  criterion difference using the restricted-fit spectrum.
* `wald_test(..., vcov = "sandwich")` constructs a complete sandwich for
  ML/FIML and continuous or all-ordered LS. Estimated weights and their joint
  first-stage influence are included; expected curvature and fixed weights are
  explicit comparators. `vcov = "fitted"` retains the original ML/FIML adapter.
  Ordinal inference requires complete data and nonempty pairwise cells. Its
  numerical reference path differentiates the actual two-step estimator and
  lavaan's native weight functional, with refined polychoric roots.
* `profile_score_test_matrices()` exposes nuisance projection separately from
  efficient-score quadratic evaluation.

* `score_test_matrices()` and `wald_test_matrices()` construct inspectable
  `semTests_test` results from supplied efficient-score matrices or parameter
  estimates and sampling covariance. `pvalues()` now dispatches on these results
  as well as lavaan fits, preserving existing lavaan defaults and output names.
  Matrix score tests support mixture calibrations; fully studentized score and
  Wald tests use their chi-square reference. These functions do not extract
  matrices from fitted models or establish misspecification validity.

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
