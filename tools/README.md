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
- classical nested method 2000 with the delta restriction map, biased gamma,
  and both ML and RLS base statistics;
- continuous GLS and ULS single-model spectra, with lavaan's explicit
  `N-1` covariance scaling and estimator-specific empirical/normal-theory
  gamma convention;
- observed-information FIML, single-model and nested (`delta` and `exact`);
- all-ordinal DWLS, ULS, and full WLS, single-model and nested, with complete and
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
- magmaan does not yet implement semTests' nested Du--Bentler `UG` spectrum;
- magmaan does not yet expose the complete-data method-2001 difference
  spectrum;
- semTests' FIML parity target is `fiml.convention = "observed"`, not its
  lavaan-compatibility convention;
- continuous GLS/ULS is certified only for single-model tests; magmaan's
  continuous-LS profile LRT is not the Satorra-2000 restriction-map spectrum
  needed for nested semTests parity.

High-throughput magmaan simulations can replace semTests only inside this
boundary. Treat an omitted combination as unverified, not as implicitly
equivalent.

The script is under `tools/`, rather than `tests/testthat/`, because magmaan is
not a CRAN package and must not become an automatic package-check dependency.
Override tolerances only for diagnosis. The optimizer-endpoint tolerances are
deliberately separate from the strict fixed-input and spectrum tolerances:

```sh
TRANSFORM_TOLERANCE=2e-6 SPECTRUM_TOLERANCE=1e-4 \
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
