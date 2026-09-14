# Non-CRAN validation

These developer checks complement the fast `testthat` suite. Run them from the
package root after installing the required local packages.

## Independent magmaan parity

```sh
Rscript tools/magmaan-validation.R
```

This opt-in check loads the current semTests source tree and compares it with
magmaan's independent C++23 implementation. It covers:

- every shared p-value transform (`STD`, `SB`, `SS`, `SF`, `ALL`, `PALL`,
  `EBA`, `pEBA`, and several `pOLS` penalties) on identical fixed spectra;
- classical complete-data ML, including biased and Du--Bentler unbiased gamma,
  ML and RLS base statistics, and the full shared test family;
- classical nested method 2000 with the delta restriction map, biased and
  Du--Bentler unbiased gamma, both ML and RLS base statistics, and a
  multigroup mean-structure check;
- continuous GLS and ULS single-model and nested spectra, with lavaan's explicit
  `N-1` covariance scaling and the estimator-specific empirical/normal-theory
  gamma convention;
- observed-information FIML, single-model and nested (`delta` and `exact`),
  for one or several groups;
- joint random-x models for continuous ML, GLS, ULS, and FIML, plus mixed
  categorical DWLS, in single-model and nested form;
- all-ordinal DWLS and ULS, single-model and nested, with complete and
  pairwise-missing data;
- mixed continuous/ordinal DWLS, single-model and nested;
- multigroup all-ordinal DWLS, single-model and nested.

The fixed-spectrum layer isolates formula parity from model fitting and uses a
strict tolerance. Spectrum and base-statistic checks then validate the
estimator-specific inputs. End-to-end p-values use separate tolerances because
lavaan and magmaan optimize the same model independently; ordinal optimizer
endpoints differ more than the spectra do, especially in multigroup fits. For
pairwise data, magmaan uses `pd_gamma = "overlap"`, matching lavaan 0.7-2's
observation-overlap convention.

The script prints an MD5 fingerprint of the installed magmaan R database and
shared library. This identifies the exact installed code even though magmaan's
development version is still `0.0.1`. Set `PARITY_OUTPUT` to retain a
machine-readable record:

```sh
PARITY_OUTPUT=workspace/magmaan-parity.csv \
  Rscript tools/magmaan-validation.R
```

### Certified simulation boundary

Passing the script certifies only the combinations above. In particular:

- complete-data nested parity requires `method = "2000"` and magmaan
  `A.method = "delta"`; magmaan's `"exact"` default is a different restriction
  construction;
- magmaan does not yet expose the complete-data method-2001 difference
  spectrum;
- semTests' FIML parity target is `fiml.convention = "observed"`, not its
  lavaan-compatibility convention.
- fixed or conditional observed exogenous predictors are outside the certified
  boundary; the random-x checks use `fixed.x = FALSE` and
  `conditional.x = FALSE`.

High-throughput magmaan simulations can replace semTests only inside this
boundary. Treat an omitted combination as unverified, not as implicitly
equivalent.

The script is under `tools/`, rather than `tests/testthat/`, because magmaan is
not a CRAN package and must not become an automatic package-check dependency.
Override tolerances only for diagnosis. The optimizer-endpoint tolerances are
deliberately separate from the strict fixed-input and spectrum tolerances:

```sh
TRANSFORM_TOLERANCE=2e-6 SPECTRUM_TOLERANCE=1e-4 \
  MULTIGROUP_SPECTRUM_TOLERANCE=5e-4 \
  PVALUE_TOLERANCE=1e-3 ENDPOINT_PVALUE_TOLERANCE=2e-3 \
  MULTIGROUP_PVALUE_TOLERANCE=5e-3 \
  Rscript tools/magmaan-validation.R
```

## Seeded null simulations

```sh
CORE_REPS=500 STRESS_REPS=300 N=500 NCORES=8 \
  Rscript tools/categorical-validation-simulation.R

NREP=500 N=400 NCORES=8 \
  Rscript tools/fiml-validation-simulation.R
```

These scripts estimate rejection rates at 1%, 5%, and 10% under seeded null
scenarios. They write CSV results under `workspace/` and report the number of
finite replications for every method. Use them to catch gross calibration
errors that deterministic equality tests cannot reveal. Small smoke runs check
execution only; hundreds of replications are needed before interpreting the
rejection rates.

## Inference implementation and timing

The phased [inference plan](score-wald-plan.md) covers computation reuse,
reduced matrix algorithms, and score/Wald extensions.

Run `Rscript tools/inference-timing.R` to compare call-local sharing with the
pre-sharing Git revision; `REPS` controls repetitions. It checks numerical and
metadata agreement before timing and writes results under `workspace/`.
Run `Rscript tools/magmaan-score-wald-validation.R` for optional score/Wald
primitive checks. These check formulas using magmaan ingredients, not yet an
independent lavaan adapter or simulation calibration.
The script also checks the exported `score_test_matrices()` and
`wald_test_matrices()` interfaces and their `pvalues()` methods against magmaan.

Run `OPENBLAS_NUM_THREADS=1 Rscript tools/residual-spectrum-timing.R` to compare
full and reduced single-model spectra on identical fitted models. `P` sets
comma-separated indicator counts divisible by three; `REPS` sets repetitions.
This excludes fitting time and saves numerical differences, timings and session
information under `workspace/`. Small models may be slower on the reduced path.

Run `OPENBLAS_NUM_THREADS=1 Rscript tools/projected-gamma-timing.R` to compare
tiled raw moment projection with uncached Gamma construction and cached Gamma
reuse, for three-dimensional and full residual-space projections. Fitting and
eigendecomposition are excluded. `P` and `REPS` work as above; results and session
information are saved under `workspace/`.

`Rscript tools/magmaan-lavaan-inference-validation.R` compares the complete-data
ML `score_test()` adapter against independent magmaan fits in five configurations.
It also checks `wald_test()` against magmaan's primitive using the same fitted
lavaan covariance. Results and installed-code fingerprints are saved under
`workspace/`; this is a numerical parity check, not a calibration simulation.

## Complete criterion inference

```sh
Rscript tools/magmaan-fiml-score-validation.R
Rscript tools/magmaan-criterion-validation.R
INFERENCE_REPS=200 Rscript tools/inference-calibration.R
```

The FIML runner checks independently extracted global scores for expected and
observed saturated projection, single and unequal groups. The criterion runner
checks continuous GLS/WLS/ULS complete and fixed-weight covariances, binary DWLS
complete Wald and restriction-release score tests, and a three-category case.
Continuous covariance parameters are explicitly converted between magmaan's N
and lavaan's N-1 sample-covariance conventions. The LS quadratic uses lavaan's
N-1 criterion scale; full studentization cancels that scalar in one group.

**First-stage convention:** magmaan's installed ordinal object stores an OPG
first-stage influence. The actual polychoric estimator has derivative
sensitivity; semTests differentiates that estimator. The three-category fixture
therefore records the native covariance gap and then supplies the derivative
influence to magmaan to isolate its independently implemented weight correction.
This is an aligned-input comparison for that fixture, not independent extraction
parity. The fitted moments and weight map agree. Binary fixtures use the native
magmaan influence without replacement. No magmaan source or installed code is
modified. CSVs distinguish these checks and fingerprint the installed binary.

The calibration runner uses an analytic misspecified-mean ML pseudo-slope of
0.4, and two independent identical ordinal populations with a misspecified
one-factor DWLS model. Equality of the pseudo-loadings follows from identical
populations and the shared criterion. It records every attempt and failure,
reports Monte Carlo standard errors, and compares sandwich Wald/score,
curvature-weighted score and criterion difference. This small preflight does
not certify a paper's design, pEBA calibration, or performance across settings.

`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript tools/criterion-inference-timing.R`
compares casewise recomputation with cached pattern-sufficient-statistic gradients
for the observed ML/FIML Hessian on identical fits. It validates the Hessian
before timing and excludes data preparation, fitting and tail evaluation.
