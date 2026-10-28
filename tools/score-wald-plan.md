# Efficient inference, score and Wald implementation plan

## Status update, 2026-10-28

B4 is settled: keep the unstructured saturated information in the released
FIML path; see the result under B4 below. The public `weight` argument of
`semtest()` keeps its name, and the two unsupported metric combinations
(observed weight with expected-curvature profiling; different score
covariances in the statistic and the reference) stay out.

**Nested restriction, decided: lavaan's `"delta"` stays the default.** The
same-point construction made the default on 2026-10-01 (e07170e) turned out
not to be invariant to identification: on the same data, `std.lv = TRUE` and
marker loadings gave nested p-values about 2e-3 apart at n = 300 (ML, MLR
with means, DWLS and FIML alike), while lavaan's mixed-point `"delta"` agreed
to 1e-6. The same-point map is exact only in one parameterization while the
information and score covariance sit at the H1 estimate, so it trades one
mixed point for another at the same `O(n^-1/2)` order. It was removed, with the
`"legacy"` name. `A.method` again takes lavaan's names: `"delta"` (default,
as in 1.0.0) and `"exact"` (key-aligned parameter bases, which shares the
identification dependence). The one kept fix from e07170e: `A.method` now
applies to complete-data fits, so `"exact"` works there too. magmaan's own
default `"exact"` has the same dependence.

The magmaan runners pair the constructions by name again (semTests `"delta"`
with magmaan `"delta"`, `"exact"` with `"exact"`). After the change (1f0fc0b),
against magmaan main b644dc56 (its pairwise `pd_gamma` fix), every runner
passes: the gate's 97 checks (824 values), the score, lavaan-adapter and
FIML-score runners, and all 20 criterion checks. The 200-draw parity
simulation (N = 500) passes every draw of all 8 cells for `"delta"`, at the
1.0.0 optimizer-endpoint floor (1.56e-3 single, 7.4e-4 nested; multigroup
ordinal 4.48e-3 and 2.17e-3), including the pairwise ordinal cell. For
`"exact"`, complete-data GLS and ULS match at the same floor; FIML differs by
up to 8.9e-3 (66.5% of draws within 2e-3) because lavaan is fitted with
`std.lv = TRUE` and magmaan with marker loadings, and `"exact"` depends on
that; categorical pairs do not offer `"exact"`. Two runner changes were needed. The
criterion runner's continuous WLS case lets magmaan build its own ADF weight,
since magmaan refuses a caller-supplied `W` for estimated-weight inference
(agreement 1.3e-7). And its ordinal nested score checks compare each
curvature with its counterpart: magmaan's ordinal LS score tests default to
observed sensitivity (`bread = "observed"`, magmaan 189fa24a), and magmaan's
`mi` and scaling factor weight by expected information whatever the bread,
which in semTests is `metric = "expected"`. Matching settings agree to about
1e-8.

## Status update, 2026-09-29

The Wald test has been removed from the package: `wald_test()`,
`wald_test_matrices()`, their tests, and the Wald checks in the `tools/`
runners. Score and criterion tests are distinct tests with their own options
and stay behind `semtest()` unchanged. The complete parameter sandwich that
the removed Wald adapter exposed is still built by `parameter_inference_parts()`
for the nested tests, and `tools/magmaan-criterion-validation.R` still checks it.
The vignette is now `vignettes/score-criterion.Rmd` and is part of the build.

B3 below was decided differently: `fiml.convention` is kept for one release
with a deprecation warning (class `semTests_warning_deprecated`) rather than
removed, because the released 1.0.0 NEWS pointed users at
`fiml.convention = "lavaan"`.

The compiled magmaan R package was renamed `magmaanlab` (magmaan commit
1866196e), and its estimate-only `magmaan()` is now `fit_model()`. The runners
now call `magmaanlab`; with that rename alone the 97-check / 824-value gate,
the score, lavaan-adapter, FIML-score and criterion runners, and a 5-draw
parity-simulation smoke run all pass. The three-category ordinal covariance
now agrees without the aligned first-stage substitution (gap 3.2e-8), because
the closed-form influence uses the same outer-product bread as magmaan.

References to Wald below are historical.

## Current implementation status

Interface naming: `curvature` selects the sensitivity matrix and
`score_covariance` selects the score covariance, with levels `"observed"` /
`"expected"` and `"centered"` / `"uncentered"`. These replace the earlier
development names (`sensitivity`/`covariance`, then `bread`/`meat`, and the
`"opg"` level) throughout the inference APIs, matrix inputs, and
ingredient/provenance fields. Names describe the matrix rather than borrowing
sandwich-literature jargon, because the audience reads lavaan documentation.
Defaults and numerics are unchanged by the rename. Wald sampling covariance
remains `vcov`; score weighting remains `metric`.

Open naming work, tracked but not yet applied: `metric` levels still read
`"information"` / `"expected"` / `"sandwich"`, which do not say which matrix
they select and collide across functions (`"information"` resolves to expected
Fisher in `score_test()` but to the selected `curvature` in
`score_test_nested()`). `fiml.convention` is slated for removal. See the
handoff notes at the end of this document.

The three remaining implementation workstreams are now implemented in the
explicitly documented scope:

* ML/FIML global score extraction, expected/observed saturated sensitivity and
  metric, and joint nested affine score releases at a restricted fit.
* Full observed-criterion curvature and complete Gaussian score covariance,
  including misspecified means, for pseudo-parameter Wald and nested score
  inference. Actual nested criterion differences share the restricted-fit
  profile spectrum.
* Joint moment/estimated-weight influence for continuous GLS/WLS/DWLS/ULS and
  complete all-ordered DWLS/WLS/ULS. Fixed weights and expected curvature remain
  explicit comparator choices.

The package remains pure R. New public entry points are
`score_test_nested()`, `criterion_test_nested()` and
`profile_score_test_matrices()`; `score_test()` now includes FIML, and
`wald_test(..., vcov = "sandwich")` constructs a new covariance. Positional
restrictions, normalization and inspectable results retain the matrix API's
contract. The executable `vignettes/score-wald.Rmd` explains use and targets.

Remaining extensions are mixed-indicator/pairwise-missing first stages,
sampling weights/clusters, nonlinear or differently parameterized nesting,
and faster ordinal first-stage differentiation. They are excluded by explicit
checks. Existing fit-test support is preserved. Raw projected Gamma for the
older FIML/categorical fit-test paths remains an optional optimization.
Paper-scale simulations and claims about penalized calibration remain owned
by the papers; the library includes only a calibration preflight.

### Algorithm and convention checks

`R/inference_geometry.R` evaluates Gaussian scores by observed pattern and
projects before covariance accumulation. Observed parameter curvature uses
central differences of aggregate analytic scores, reusing pattern sufficient
statistics at all perturbations. This avoids repeated casewise meat extraction.
`R/inference_weights.R` differentiates continuous Gamma analytically, including
sample-mean effects. GLS and continuous DWLS differentiate their native
normal-theory sample-covariance weight map. The LS score, curvature and influence
all retain lavaan's per-group N-1 criterion factors.

The ordinal reference path differentiates the native frequency-weighted
functional, using locally refined polychoric score roots and Richardson
extrapolation. It does not modify lavaan's namespace or refit the structural
model. Duplicate response patterns share a derivative; DWLS retains diagonal
weight derivatives rather than full matrices. Empty pairwise cells are rejected.

The multi-category check exposed a first-stage convention difference in
magmaan: its stored influence uses OPG bread, while the derivative of the actual
polychoric estimator uses derivative sensitivity. Both packages retain the same
OPG-based weight map. In the three-category fixture, the native covariance gap
is about .00141. Supplying the actual estimator influence to magmaan isolates
its independent weight calculus and reduces the covariance gap below 3.1e-8.
The binary fixture agrees independently without substitution (below 5e-9).
This distinction is retained in the runner and diagnostics CSV, rather than
hidden by a loose tolerance. No magmaan files were changed.

`tools/magmaan-criterion-validation.R` passes 19 checks, covering continuous
complete/fixed covariance and Wald, plus binary and aligned-input multi-category
DWLS covariance and nested scores. The full continuous covariance gaps are
below 1.3e-7 after the explicit N versus N-1 parameter conversion.
`tools/magmaan-fiml-score-validation.R` passes 36 checks of saturated sensitivity and
metric choices, including matching rejection of an indefinite metric.
The original 97-check / 824-value magmaan gate remains green.

The unit suite checks likelihood and Hessian derivatives, coordinate invariance,
FIML reduction to complete ML, pattern tiling, affine template embedding,
unequal-group criterion normalization, empirical contamination derivatives,
first-stage step convergence and analytic threshold influences. Validation
reports are described in `tools/README.md`. The focused criterion file now
passes 69 assertions. The final source build and `R CMD check --no-manual`, including installed
tests, examples and all seven vignettes, passed with Status: OK after the
affine starting-value fix. The final log is
`/tmp/semtests-complete-inference-check/final-check.log`.

The template adapter uses the alternative's actual affine constraint offset,
not its potentially infeasible starting estimates. ML/MLM/MLR labels may differ
between restricted fit and template when they define the same criterion.
A three-repetition, single-BLAS-thread timing check measured 0.301 to 0.030 s
(complete ML) and 0.338 to 0.032 s (FIML) for the 30-parameter Hessian. Differences
were below 1.5e-9. This is curvature-only timing, not end-to-end inference speed.

### Calibration preflight

`tools/inference-calibration.R` completed 200 attempts in each of two designs,
with zero convergence/admissibility/curvature failures. It retains every
attempt in `workspace/inference-calibration-attempts.csv` and reports uncertainty
in `workspace/inference-calibration-summary.csv`. Seed: 420190. These are small
library preflights, not a frozen paper simulation bank.

| Method | ML misspecified means, N=600 | Ordinal DWLS, N=400 per group |
| --- | --- | --- |
| Complete Wald | 6.5% (MCSE 1.74 pp) | 8.5% (MCSE 1.97 pp) |
| Sandwich score | 6.0% (MCSE 1.68 pp) | 2.5% (MCSE 1.10 pp) |
| Curvature-weighted score | 6.0% (MCSE 1.68 pp) | 3.5% (MCSE 1.30 pp) |
| Criterion difference | 6.0% (MCSE 1.68 pp) | 3.5% (MCSE 1.30 pp) |

All rejection rates use nominal 5%. The ML pseudo-slope is analytically .4
because E[XY]/E[X^2]=.8/2, while both fitted means are incorrectly fixed at
zero. The ordinal null follows from two identical independent populations
sharing the DWLS criterion, with four equal pseudo-loadings and a misspecified
one-factor model. The ordinal point estimates show conservative score and
liberal Wald behavior in this finite sample. These results do not establish
nominal size, superiority, or pEBA calibration; larger paper-owned experiments
remain necessary for such claims.

The following execution plan and phase descriptions preserve the incremental
implementation history; this current-status section supersedes their earlier
statements about work still to be implemented.

## Execution plan

Keep the implementation in R initially. Carry over the algebraic reductions
from magmaan, independently verifying the fitted-model conventions. Add C++
only for a measured bottleneck remaining after these reductions.

| Phase | Deliverable | Acceptance gate |
| --- | --- | --- |
| 1: shared computation | Call-local reuse of statistics/spectra, standard-only bypass, reproducible timing harness | Batch/single agreement, ordering and suffix checks, package tests and existing magmaan parity |
| 2: smaller matrix problems | Shared symmetric factor-solve kernel; expected-information single-model residual-space spectrum; retain current route as reference | Fixed-input and fitted-model parity, unequal groups and constrained models, scaling benchmarks |
| 3: projected covariance | Blockwise/tiled moment influence projection, avoiding full Gamma where valid | Compare projected products with materialized Gamma; separate missing/categorical/unbiased coverage |
| 4: score/Wald foundation | Matrix kernels, common inspectable result, then ML and FIML score adapters and linear Wald | Independent derivatives and magmaan primitives, no-missingness reduction and restriction invariance |
| 5: fixed misspecification | Full observed curvature and complete NTML score covariance | Resolve mean-score omission, independently established pseudo-null tests |
| 6: estimated weights | DWLS joint first-stage influences, then GLS/WLS | Full/fixed-weight comparators and exact-fit reduction; paper-owned simulation designs |

Each phase is a reviewable increment. Do not silently extend a fast path to an
estimator or convention that has not passed its gate. Keep extraction, linear
algebra and tail evaluation timings separate. Retain backwards compatibility
for the existing public calls throughout.

Phase 1 is implemented in `R/pvalues.R`. The cache is scoped to a
single public call, so changing fitted objects or conventions cannot reuse stale
ingredients. Spectrum reuse is independent of the requested base statistic.
`tools/inference-timing.R` compares against the pre-sharing source revision and
checks identical results including metadata. A 10-repetition HS CFA run measured
1.86x for SB/SS/SF with ML, 2.06x with RLS, and 1.02x for SB/ALL/pEBA4 with ML.
These local timings support sharing but also show that mixture-tail work needs
separate attention. The existing magmaan gate passed 97 checks / 824 values.

Phase 2 is implemented in `R/quadratic.R`, `R/gamma.R` and `R/fiml_fmg.R`:
expected-information single-model residual factors, groupwise Gamma products,
the observed-convention FIML residual factor when its saturated metric is PD,
and symmetric Cholesky solves for nested spectra. Observed model information
outside this FIML construction retains the lavaan route. Non-PD metrics retain
the original algorithm; rank failures in the reduced calculation are explicit.
That phase still extracted full Gamma; phase 3 adds direct influence projection.

`tools/residual-spectrum-timing.R` checks the same fitted objects against full
lavaan UGamma before timing. Initial single-thread runs (three repetitions)
gave 1.22x at 18 indicators and 1.58x at 30, with maximum eigenvalue differences
below 7e-14. The nine-indicator example was slower (3 ms versus 1.7 ms), so this
is a scaling improvement, not a claim that every model is faster. The unchanged
magmaan parity gate passed all 97 checks / 824 values after this change.

Phase 3's complete-data ML slice is implemented in `R/moment_influence.R`.
`projected_gamma()` reuses cached NACOV or accumulates projected empirical
moment influences in row tiles when NACOV is absent. It supplies both single-
model and nested biased spectra. All biased non-FIML nested paths now accumulate
group covariance products without a global block-diagonal Gamma matrix.
The raw route requires complete-data normal-likelihood ML with the standard
Gamma normalization, no sampling weights/clusters, and no supplied NACOV or
sample covariance. Means and unequal group sizes are included. Unbiased Gamma,
other estimators and alternative first-stage conventions retain their verified
materialized covariance routes. FIML and categorical raw influences remain
separate future work; phase 3 does not claim to implement those conventions.

The raw kernel uses O(n*p) storage for centered observations plus row-tile
moment products and the projected covariance. It constructs neither a full
n-by-m moment-influence matrix nor an m-by-m Gamma. This reduces additional
post-fit memory; it cannot remove Gamma that lavaan already built during fitting.
The cached path is intentional: recomputing raw influences can cost more than
multiplying a covariance that is already available. Forced raw/materialized
sources are internal validation and timing controls, not public API options.
The independent tiled tests cover mean/covariance cross terms, scalar and
multidimensional factors, uneven tiles, translation invariance, unequal groups,
uncached single/nested fits, and custom NACOV preservation.

Validation: 34 focused assertions and the existing 97-check / 824-value magmaan
gate pass. Projection-only timings (`tools/projected-gamma-timing.R`, 500 cases,
one BLAS thread, five repetitions) at 30 indicators were 12.4 ms for uncached
Gamma versus 2.8 ms for tiled projection into three dimensions (4.4x). For
the full 402-dimensional residual space they were 30.4 versus 16.6 ms (1.8x).
These exclude fitting and eigendecomposition and are not end-to-end speedups.

Phase 4's matrix foundation is implemented in `R/inference_tests.R`:
`score_test_matrices()` consumes already efficient scores, metric and covariance;
`wald_test_matrices()` consumes theta, a linear restriction and theta's sampling
covariance. Both return a `semTests_test` with statistic, rank/df, spectrum,
reference, hypothesis, normalization, diagnostics, provenance and optional
ingredients. Print and data-frame methods support inspection and saved results.
`pvalues()` is now an S3 generic, with the original lavaan implementation retained
as its lavaan method. Existing lavaan defaults and output names are preserved.

Matrix score inputs explicitly distinguish summed versus mean/per-observation
units. The supplied-metric default is the full mixture (`ALL`); setting the
metric to NULL fully studentizes with the supplied covariance. Fully studentized
score and Wald results use `STD` only. This prevents a second eigenvalue
correction. Matrix results reject lavaan statistic/Gamma suffixes. Wald never
adds a second N factor to an already supplied sampling covariance.

The initial matrix interface rejects singular score metrics/covariances and
redundant or unidentified Wald restrictions. Theta's full covariance may be PSD
when the tested restriction covariance remains PD. The package does not infer
the correctness of supplied derivatives or a misspecification target.
The 48-check magmaan primitive script now exercises the public APIs as well as
independent formulas; the existing 97-check parity script still passes.
The new unit file passes 49 assertions. A source build with all six vignettes
and `R CMD check --no-manual` completed with `Status: OK`, including installed
tests, examples, documentation, method registration and vignette rebuilds.
The explicit `_matrices` names distinguish this foundation from fitted-model
adapters and keep independently constructed covariance estimators usable.

Phase 4 now also includes complete-data ML adapters in `R/inference_lavaan.R`.
`score_test()` builds direct saturated Gaussian log-likelihood scores at the
fitted mean/covariance and removes the model tangent using expected Fisher
information. It tiles observations, returns a Fisher-whitened residual basis,
and supports OPG or group-centered covariance and optional full studentization.
It includes mean/covariance score cross terms, unequal groups, affine equalities
and joint random-x models. The null is global exact fit, not equal population
fit under fixed misspecification. Saturated observed-sensitivity variants and
nested score releases remain unimplemented.

`wald_test()` uses `lavaan::coef()` and the fitted sampling covariance, recording
parameter order and the fitted SE/information settings. It rejects restrictions
already fixed or dependent in the model's effective parameter space. This is
an adapter to the fitted covariance, not a newly constructed complete sandwich.

Validation adds 46 unit assertions, including finite-difference likelihood
scores, a population likelihood Hessian check, tangent orthogonality, tiling,
and agreement with `lavaan::lavTestWald()` for ML/MLM/MLR. The new optional
`tools/magmaan-lavaan-inference-validation.R` passes 45 checks across covariance-
only, mean-structure, constrained-mean, unequal-group/equal-loading and random-x
fits. Its score tests use independently fitted and extracted lavaan/magmaan
objects; Wald compares with the C++ primitive on the same fitted covariance.
The existing 97-check / 824-value parity gate also passes. A full source build
and `R CMD check --no-manual`, including installed tests, new examples and all
vignettes, completed with `Status: OK` for these adapters.
FIML extraction and nuisance projection are the next part of phase 4, followed
by the distinct fixed-misspecification and estimated-weight workstreams.

Initial planning snapshot: 2026-09-12. The phase status above records subsequent
implementation. Initial semTests source HEAD: e54d90b27c03f151a73ad57a456308ad44862328;
magmaan source HEAD: 13c7e69c28dcd7c2efb727f216b0452d21123b8a.
The installed magmaan binary is fingerprinted by the existing parity runner;
its identity must not be inferred from the source checkout's HEAD.

## Located implementations and paper requirements

Reference checkout: `/home/jonas/Files/research/magmaan`.
The papers are `papers/fiml-fmg` and `papers/sem-misspecification`.
Read their current AGENTS.md before relying on older manuscripts or simulations.
The FIML paper now concerns score/LRT/pEBA; its legacy invariance design is stale.
The misspecification paper includes NTML and categorical DWLS, pseudo-parameter
restrictions and equal population fit. Missing-data extensions are optional there.

Verified source surfaces:

| Surface in magmaan | What it establishes |
| --- | --- |
| `global_score_flip_test()`; `include/magmaan/inference/score.hpp` | Global ML/FIML/ML2S score construction, asymptotic mixture and sandwich comparators, projected score/metric/meat. Observed sensitivity and metric options have narrower support. |
| `score_tests_robust()` | Candidate-wise restriction-release tests; exposes bread and estimated-weight choices. This is not automatically a joint global test. |
| `weighted_param_space_sandwich_ij()` in `robust/weighted_inference.hpp` | Estimated-weight contribution to score covariance. Its documented bread is the expected moment-Jacobian bread; a complete observed-curvature joint score path needs separate verification. |
| `infer_ordinal_robust_ij()` and mixed counterpart | Complete estimated-weight sandwich ingredients exist for parameter inference. |
| `infer_wald_test_theta(theta, R, vcov, q)` | General linear Wald quadratic using a supplied sampling covariance. It does not construct that covariance or a global-fit restriction. |
| `*_profile_lrt()` | Separate criterion-difference reference constructions; not interchangeable with score or Wald statistics. |

Thus estimated-weight Wald is possible through a complete covariance plus the
Wald primitive. A dedicated general weighted global Wald API was not established
by this inspection. Eigenvalue weights, estimation weights, and observation
weights must be named separately. This plan covers the first two; survey/case
weights need their own specification.

## Statistical contracts

1. **Exact fit:** test the fitted model against saturated moments. Distributional
   robustness allows nonnormality/missingness under the stated assumptions. An
   observed-information option alone does not establish fixed-misspecification validity.
2. **Pseudo-parameter restrictions:** under fixed misspecification test
   `r(theta*) = 0`, with theta* defined by the actual population criterion.
   Use its full observed curvature and full estimating-function covariance.
3. **Equal population criterion:** nested fits must share the criterion and
   estimated-weight map. Implement the criterion-difference expansion separately.

Global exact fit is false under fixed structural misspecification. Calling a
test "misspecification robust goodness of fit" must not imply nominal null size
there. A global Wald test additionally needs an explicit saturated-coordinate
restriction map; arbitrary parameter Wald tests do not provide one.

For efficient score vector s in summed units, metric H and score covariance J
in matching summed units, the common quadratic is s' H^-1 s and its reference
spectrum is eig(H^-1/2 J H^-1/2). The fully studentized score uses s' J^-1 s.
For Wald, d = r(theta_hat) and its sampling covariance V_r give d' V_r^-1 d;
do not multiply by N again when V_r is already a sampling covariance.
Full sandwich studentization and an arbitrary metric with a mixture reference
are different choices. Do not apply pEBA a second time to an already whitened test.

Estimated discrepancies require joint first-stage moment/weight influences,
including their cross covariance and derivative terms. Retain residual-curvature
terms in the bread. Fixed-weight and expected-bread variants are useful explicit
comparators, not silent fallbacks for the complete correction.

## Recommended library shape

Keep semTests an independent R package over lavaan, with magmaan an optional
development oracle. Avoid adding a compiler requirement or runtime magmaan
dependency before profiling identifies a bottleneck.

Algorithmic improvements are in scope regardless of implementation language.
The current package has no `src/`, Rcpp import, or LinkingTo dependency. R's
matrix products and decompositions already use compiled BLAS/LAPACK.

Optimize the existing tests alongside the new kernels:

* Implemented in phase 1: extract statistics and spectrum once per distinct
  fit/convention/gamma choice in a call, then apply all requested calibrations.
  Standard-only requests skip spectrum construction. `compute_pvalues()` still
  handles each parsed test but shares expensive ingredients through the context.
* Keep the existing restriction-space reduction for nested tests. Replace
  nonsymmetric inverse-products with symmetric whitened eigenproblems where
  definiteness permits, using factor solves and explicit rank checks.
* Add a single-model factor path: where the actual construction admits
  `U = B B'`, diagonalize `B' Gamma B` in residual dimension instead of the
  full moment-space `U Gamma`. magmaan's `finish_u_factor_expected()` constructs
  the complement using QR without materializing the full orthogonal matrix.
  Preserve unequal group weights and equality constraints. Do not assume this
  PSD factorization for every observed-bread or misspecification construction.
* For empirical `Gamma = Z'Z/N`, project moment-influence rows first and use
  `crossprod(Z %*% B)/N`, accumulating tiles when needed. Avoid constructing
  full Gamma when only its projected covariance is needed. Missingness,
  categorical first stages and unbiased Gamma need their own equivalent paths.
* Avoid dense block-diagonal group matrices and explicit Kronecker products;
  accumulate blockwise products and use covariance-index formulas directly.
* Profile extraction, statistic construction (including RLS), projection,
  eigendecomposition and tail evaluation separately. Optimize the R Ruben
  recurrence's repeated vector growth before deciding whether a small compiled
  tail kernel is justified. FIML pattern-score loops are another possible
  compiled target, not a reason to migrate the whole package.

A small local timing check (100 repeated calls, nine-variable HS CFA, MLM,
OpenBLAS) gave approximately 0.97 ms for `lavInspect(fit, "UGamma")`, 0.31 ms
for its 45-dimensional eigenspectrum, and 192 ms for a public call requesting
SB/ALL/pEBA4 with default statistic selection. These are illustrative timings,
not a scaling benchmark or an attribution of the remaining time to one stage.
Benchmark larger models and isolate tail/statistic costs before selecting C++.

Introduce three layers:

* **Extraction:** `R/inference_lavaan.R` and `R/inference_fiml.R` produce a validated
  context: parameter/moment ordering, restriction/tangent bases, residuals,
  Hessian, casewise or first-stage influences, group sizes and conventions.
  Reuse the alignment and comparability checks already in `fiml_fmg.R` and
  `lavaan_helper.R`; do not reuse an LR spectrum as a score spectrum.
* **Numerical kernels:** `R/score.R`, `R/wald.R`, `R/quadratic.R`, and later
  `R/estimated_weight.R`. Matrix inputs permit independent unit tests and external
  adapters. No fitting, estimator dispatch, or package-specific objects here.
* **Public results:** proposed `score_test()` and `wald_test()` return a common
  `semTests_test` object; `pvalues()` gains S3 dispatch while preserving existing
  lavaan calls, defaults, return names and metadata. This requires an explicit
  compatibility test before conversion of the current ordinary function.

Illustrative API, subject to the first implementation slice:

```r
s <- score_test(fit, hypothesis = "overall", bread = "expected",
                metric = "information", meat = "empirical")
pvalues(s, tests = c("SB", "ALL", "pEBA4"))
w <- wald_test(fit, R = R, rhs = q, bread = "observed",
               meat = "empirical", weight_effect = "estimated")
as.data.frame(w)
```

Keep sensitivity, test metric, distributional covariance and estimated-weight
effect separate. No omnibus `robust=TRUE` switch. Arbitrary-metric Wald and
global Wald should remain deferred until their restrictions and target laws are
specified. Initial Wald supports linear restrictions and full sandwich inference.
Initial score needs only the restricted fit plus the alternative specification
or restriction map, not an unnecessary fitted alternative.

Result fields: statistic, nominal df, effective rank, eigenvalues, reference law,
calibration results, hypothesis/target, estimator, sensitivity, metric,
covariance/centering convention, weight treatment, normalization and sample sizes.
Include convergence/stationarity, minimum eigenvalues, condition diagnostics,
source versions and numerical tolerances. Retain ingredients optionally for
research inspection. Fail explicitly on unsupported combinations, material
negative eigenvalues or unidentified restrictions; no silent pseudoinversion
or truncation to nominal df. Signed criterion spectra need a separate law path.

## Implementation order and acceptance gates

1. **Matrix kernels and result contract.** Score quadratic/spectrum and linear
   Wald; fixed-input agreement with magmaan; independent transformations,
   scaling and restriction-basis invariance tests. Separate H, J and V naming.
2. **Complete ML exact-fit global score.** Independent lavaan extraction, same
   model/data against magmaan, compare score ingredients and spectra before
   p-values. Add nested linear score restrictions next.
3. **FIML score.** Observed-data pattern scores, explicit expected versus observed
   sensitivity/metric, group normalization and nuisance projection. Test no
   missingness reduction, MAR/MCAR examples and multigroup equality constraints.
   This is the first paper-facing milestone; reuse existing pEBA functions.
4. **Misspecified complete ML score/Wald.** Full casewise likelihood derivatives
   and observed Hessian, including misspecified means. The paper STATUS.md flags
   an installed generic NTML score helper omitting a mean-residual term. Resolve
   that by finite differences/direct likelihood scores before certifying it as
   an oracle. Establish pseudo-null targets independently of simulation draws.
5. **Estimated-weight DWLS, then continuous GLS/WLS.** Audit the exact first-stage
   weight functional and influence derivatives. Compare expected/fixed,
   observed/fixed and observed/complete corrections on identical inputs. Check
   exact-fit disappearance of residual-driven weight effects. The current blanket
   WLS refusal in semTests applies to its exact-fit correction path; reconsider
   it locally for misspecification inference rather than deleting it globally.
6. **Paper validation and documentation.** Multi-restriction tests, multigroup,
   rank-deficient cases and independently established pseudo-null simulations.
   Publish Monte Carlo uncertainty, convergence/failure denominators and pinned
   software provenance. pEBA outside its currently studied settings is an
   experimental calibration choice until evaluated, even if formula parity holds.

Do not source paper code into this library. Put portable fixtures and independent
checks here; papers consume package APIs and own their scientific designs.

## Checks executed during planning

`PARITY_OUTPUT=workspace/magmaan-parity-planning.csv Rscript tools/magmaan-validation.R`
passed all 97 checks (824 compared values). This certifies only the existing
runner's coverage; it does not certify the proposed score/Wald adapters.

`Rscript tools/magmaan-score-wald-validation.R` checks ML and FIML global score
quadratics, spectra, sandwich quadratics and mixture/SB probabilities under
expected and observed sensitivity, plus a fixed-input linear Wald quadratic.
It uses magmaan's extracted ingredients, so it tests formula agreement rather
than independent extraction or calibration. Its 19 multiplier draws exercise
the wrapper only. Outputs live in ignored `workspace/`.

No new estimator or public function is part of phase 1. Its implementation changes
computation reuse only; subsequent phases change the matrix algorithms.

## Interface handoff: remaining vocabulary work

Phase A (applied) renamed `bread` to `curvature`, `meat` to `score_covariance`,
and the `"opg"` level to `"uncentered"`, across `R/inference_*.R`, the four
inference test files, the four `tools/` validation runners, and
`vignettes/score-wald.Rmd`. Roxygen was regenerated. The suite was green before
and after at 657 passing assertions over 128 files, with the same four
CompQuadForm skips, so the rename changed no numbers and the magmaan parity
gate is unaffected. No tracked/released file was touched; every renamed file is
still untracked, so there is nothing to deprecate.

The remaining items each change numbers, so each needs the magmaan gate re-run
on its own commit. They are listed in dependency order.

### B1. Collapse the `metric` levels

`metric` currently reads `c("information", "expected", "sandwich")` in
`score_test_nested()` and `c("information", "sandwich", "observed")` in
`score_test()`. No level says which matrix it selects, and `"information"`
means expected Fisher in the first function but the selected `curvature`
(default `"observed"`) in the second. Target:

    metric = c("curvature", "score_covariance")

`"curvature"` weights by the selected curvature and keeps the weighted
chi-square mixture reference; `"score_covariance"` studentizes and takes the
chi-square reference. This drops the mixed comparator reachable today as
`metric = "expected"` with `curvature = "observed"` (expected weighting, observed
nuisance projection). Confirm nothing in `tools/magmaan-*` depends on that
combination before removing it. If an expected-weighted quadratic is still
wanted it is reachable as `metric = "curvature"` with `curvature = "expected"`.

### B2. Decide the `score_covariance` default

`score_test()` defaults to `"uncentered"`; `score_test_nested()`,
`criterion_test_nested()` and `wald_test()` default to `"centered"`. Both levels
stay for now, per the decision to preserve magmaan parity. Centered is valid in
both regimes and uncentered is consistent only when the score has mean zero, so
centered is the likely target default once the magmaan side is settled. Blocked
on confirming whether magmaan's global-score convention accumulates uncentered.

### B3. Remove `fiml.convention`

Drop the `"lavaan"` branch in `fiml_lambdas()` and `fiml_nested_ingredients()`,
which is a passthrough to `lavInspect(fit, "UGamma")` and to
`lavTech(fit, "WLS.V")` respectively. It is not a statistical convention: it
inherits lavaan's `information` and `h1.information` options, so the same call
on the same data returns different numbers depending on how the fit was
requested. Verified: `UGamma` moves with both options (5.99e-01 and 4.72e-01 on
a four-indicator FIML fixture); `WLS.V` moves only with `h1.information`, and
is unaffected by `information`, so no lavaan option yields an observed H1 weight
under FIML.

This one touches released code: `R/fiml_fmg.R`, `R/pvalues.R`, `R/support.R`,
plus four test files, `vignettes/fiml-missing-data.Rmd` and three `tools/`
runners. 103 references over 20 files, but `man/` regenerates. Approved to break
without a deprecation cycle, since the argument predates the FIML manuscript and
has no known external users. Note it in the README.

### B4. Settle structured vs unstructured saturated evaluation

`fiml_h1_information_observed()` force-overrides to `h1.information =
"unstructured"` (`R/fiml_fmg.R:214`), while the newer Gaussian inference path
builds every saturated geometry from `lavInspect(object, "implied")`, i.e. the
structured convention (`R/inference_lavaan.R:60,71`;
`R/inference_geometry.R:101,128`). Two exact-fit tests in one package currently
disagree on the same object.

Structured is defensible for a goodness-of-fit null: under the null and under
local misspecification both are consistent for the same limit, and under fixed
misspecification the statistic diverges so the reference law does not matter.
Rao's score test conventionally evaluates information at the null-restricted
estimate, which is what structured is, and the implied covariance is always
positive definite. But the naive efficiency argument for structured is not
supported by the FIML coverage table in `vignettes/fiml-missing-data.Rmd`, where
the analogous shortcut was a hair worse in its most favourable cell (94.48 vs
94.73), and that study targeted parameter coverage rather than goodness-of-fit
null calibration, so it does not settle this.

Proposed settling simulation, laptop scale: N in {100, 200, 400}; normal and
t5/Vale-Maurelli; correct model and a mild omitted cross-loading around .15-.25;
structured and unstructured; outcome is rejection rate at 5% plus power in the
transition region, for SB, pEBA4 and PALL. 24 cells. A null result is a usable
result: if the gap is inside MCSE, pick structured on principle and say so.

Note the stakes. If structured wins or ties, the path that needs changing is the
released `fiml_h1_information_observed()` override, not the new inference code,
and that is what the FIML manuscript reports.

**Result, 2026-10-28: unstructured wins; the released override stays.**
`tools/fiml-h1-information-simulation.R`, 2000 replications per data cell on
simbox (lavaan 0.7.2). MAR missingness as in
`tools/fiml-validation-simulation.R`; misspecification is an omitted .20
cross-loading of y4 on f1. Both conventions are computed on every draw, with
the same statistic, Gamma and Delta, so the comparison is paired. Switching
only the saturated information to structured makes every test reject more
often in every cell, by 4.8 to 18 paired Monte Carlo standard errors. Rejection
at 5% under the correct model, SB, unstructured vs structured:

| N | normal | t5 |
|--:|:--|:--|
| 100 | .101 vs .171 | .114 vs .256 |
| 200 | .072 vs .097 | .112 vs .175 |
| 400 | .060 vs .071 | .073 vs .105 |

pEBA4 and PALL behave the same. Size-adjusted power (critical value from the
matching correct-model cell) is equal within 1.5 points, so the extra
rejections are size distortion, not power. The simulation does not test
`semtest()`'s FIML score test, which is structured throughout (its score
covariance also sits at the implied moments), so it does not show that
`semtest()` should switch.

### B5. Open design question

Expected score covariance is wanted for FIML and goodness of fit, where it is
reported to outperform its competitors substantially. Under complete data the
information equality makes an expected score covariance the same matrix as the
expected curvature, and its goodness-of-fit consequence is already exposed as
`tests = "std"`. Under missing data the equality need not hold, so a distinct
`score_covariance = "expected"` level may be warranted. Specify which matrix is
intended before adding the level.
