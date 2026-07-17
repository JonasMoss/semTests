# Non-CRAN validation

These developer checks complement the fast `testthat` suite. Run them from the
package root after installing the required local packages.

## Independent magmaan parity

```sh
Rscript tools/magmaan-validation.R
```

This opt-in check loads the current semTests source tree and compares it with
magmaan's independent C++23 implementation. It covers:

- observed-information FIML, single-model and nested (`delta` and `exact`);
- categorical DWLS, single-model and nested;
- complete and pairwise-missing categorical data.

The strongest check is spectrum parity. End-to-end p-values use a looser
tolerance because lavaan and magmaan optimize the same model independently.
For pairwise data, magmaan uses `pd_gamma = "overlap"`, matching lavaan 0.7-2's
observation-overlap convention.

The script is under `tools/`, rather than `tests/testthat/`, because magmaan is
not a CRAN package and must not become an automatic package-check dependency.
Override tolerances only for diagnosis:

```sh
SPECTRUM_TOLERANCE=1e-5 PVALUE_TOLERANCE=5e-4 \
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
