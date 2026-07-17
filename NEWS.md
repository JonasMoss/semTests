# semTests 0.9.0

A major release on two fronts: semTests now has a self-contained numerical core
with a minimal dependency footprint (`Imports` is just `lavaan` and the base
`methods` package), and the eigenvalue-based p-values have been broadened well
beyond normal-theory ML (see `?semTests-support`).

* Adapted statistic extraction to the named-list `lavaan::lavTest()` structure
  introduced in lavaan 0.7-1 and used by the required lavaan 0.7-2 (the
  reweighted least-squares statistic was previously read as `NULL`).
* The eigenvalue-based p-values no longer depend on CompQuadForm. The upper tail
  of the quadratic form `sum lambda_j chi^2_1` is computed internally by
  `imhof_pvalue()`, which is more accurate in the tail, where CompQuadForm
  degrades (its `imhof` returns quadrature noise, including negative
  probabilities, once the p-value drops below roughly `1e-6`, and `davies`
  underflows to zero). Positive spectra use Ruben's series (exact to machine
  precision, including the tail); mixed-sign spectra use the Imhof integral in
  the body with a Lugannani-Rice saddlepoint fallback in the tail.
* Nested tests now read the reference spectrum from the reduced `m x m`
  restriction-space matrix instead of the full `q x q` UGamma, where `m` is the
  number of restrictions. The result is identical to machine precision but
  faster on large models, and it removes the need for `RSpectra`.
* Replaced the internal use of `MASS::ginv` with a self-contained
  `generalized_inverse()`.
* Removed the unused internal bootstrap helpers.
* Dropped the `CompQuadForm`, `Matrix`, `MASS`, `RSpectra`, `progressr`, and
  `future.apply` dependencies, plus the `psych` suggestion (examples now use the
  built-in `HolzingerSwineford1939` data). `Imports` is now `lavaan` and
  `methods`.
* Broadened the eigenvalue-based p-values beyond normal-theory ML. `pvalues()`
  now supports GLS, ULS, and categorical DWLS/ULS/WLS families in addition to ML/MLM/MLR,
  plus FIML missing-data fits (single-group, continuous); `pvalues_nested()`
  supports the continuous estimators, nested FIML comparison, and categorical
  Satorra-2000 comparisons with the delta restriction map. See
  `?semTests-support`. **This broadened support is experimental**; the
  classical normal-theory ML path remains stable.
* FIML now requires lavaan 0.7-2, whose corrected missing-data `Gamma` and
  observed H1 information are used directly. The former hand-coded saturated
  score and Hessian implementation has been removed.
* Added `fiml.convention` to `pvalues()` and `pvalues_nested()`. The default,
  `"observed"`, uses observed saturated and model information throughout and is
  validated against an independent magmaan implementation. `"lavaan"`
  reproduces lavaan 0.7-2's inspected single-model `UGamma` and public nested
  Satorra-2000 scaled/scaled-shifted tests.
* Biased single-model spectra now use lavaan's inspected `UGamma` directly.
  Experimental categorical support covers DWLS/WLSMV, ULS/ULSMV, and full WLS
  for ordered or mixed indicators, single or multiple groups, and listwise or
  pairwise missingness. Nested categorical tests use Satorra 2000 with the
  delta restriction map.
* Fixed constrained FIML models: equality-constraint bases are now applied to
  single-model spectra, general equality constraints mixed with inequalities
  are handled, and nested tests support an already-constrained H1.
* `A.method = "delta"` is now the nested FIML default. Exact restriction maps
  align parameters by model identity rather than parameter-table row order and
  fail clearly when the two models do not share a full parameter set.
* Removed the unreachable missing-data rescaling and the expected-information
  warning. The FIML information convention is now an explicit argument and is
  recorded in result provenance.
* Added entry-point validation with clear messages pointing at
  `?semTests-support`: non-`lavaan` objects, unsupported estimators, unsupported
  missing-data modes, multi-group / fixed-exogenous FIML, the normal-theory-only
  RLS statistic and the Du-Bentler `UG` gamma off the classical case, and
  incompatible nested pairs are now refused up front. All nested families must
  use the same estimator, fitting conventions, variables, groups, sample sizes,
  raw data, and missingness mask.
* Added typed public conditions for malformed test specifications,
  nonconverged or inadmissible fits, reversed nested inputs, incompatible model
  pairs, and unstable method-2001 spectra. Test names now enforce their grammar
  and parameter domains. Reversed inputs warn before being swapped, and method
  2001 never silently falls back to method 2000.
* The returned value is now a `semTests_pvalues` object that records the options
  actually used (requested and base estimator, requested tests, base statistic,
  information type, gamma type, data type, and degrees of freedom) and prints a
  one-line provenance footer.
* Added a `semTests` vignette and an `inst/CITATION`.
* Documentation fixes and a substantially expanded test suite.
