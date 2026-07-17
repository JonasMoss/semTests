
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semTests <img src="man/figures/logo.png" align="right" width="210" height="130" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/semTests)](https://cran.r-project.org/package=semTests)
[![R-CMD-check](https://github.com/JonasMoss/semTests/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JonasMoss/semTests/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/JonasMoss/semTests/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JonasMoss/semTests?branch=main)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

An R package for robust test statistics in structural equation models,
covering both overall goodness of fit testing and nested model
comparison. Built on top of `lavaan`. It implements the penalized
eigenvalue block averaging and penalized regression *p* values of
Foldnes, Moss, & Grønneberg (2024), together with their extension to
nested model comparison in Foldnes, Grønneberg, & Moss (2026).

The eigenvalue correction targets the limiting null law of the test
statistic, a weighted sum of chi squares. This law holds for any minimum
discrepancy estimator, so support now reaches beyond normal theory ML to
GLS, ULS, categorical DWLS/ULS/WLS families, and FIML missing data fits.
See `?semTests-support` for the full validated matrix.

## Installation

Either install from `CRAN`,

``` r
install.packages("semTests")
```

Or install the development version using the following command inside
`R`.

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/semTests")
```

## Usage

Call the `library` function, create a `lavaan` model, and run the
`pvalues` function.

``` r
library("semTests")
model <- "visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9"
object <- lavaan::cfa(model, lavaan::HolzingerSwineford1939, estimator = "MLM")
pvalues(object)
#>    peba4_rls 
#> 5.165737e-07 
#> estimator: ML (MLM) | data: continuous | information: expected | df: 24
```

For nested model comparison, fit the constrained and unconstrained
models and pass both to `pvalues_nested`, the constrained model first.

``` r
constrained <- "visual  =~ x1 + x2 + x3
                textual =~ a*x4 + a*x5 + a*x6
                speed   =~ x7 + x8 + x9"
m1 <- lavaan::cfa(model, lavaan::HolzingerSwineford1939, estimator = "MLM")
m0 <- lavaan::cfa(constrained, lavaan::HolzingerSwineford1939, estimator = "MLM")
pvalues_nested(m0, m1)
#> pall_ug_ml 
#>  0.0166158 
#> estimator: ML (MLM) | data: continuous | information: expected | df: 2 | nested (method 2000)
```

> **Note.** Support for estimators other than normal theory ML (GLS,
> ULS, categorical DWLS/ULS/WLS, and FIML) and for nested comparison
> under FIML or categorical estimation is **experimental** in this
> release. The classical ML path is stable. See `?semTests-support`.

Missing data is handled through full information maximum likelihood. Fit
the model with `missing = "fiml"` and pass it to `pvalues` as usual.
FIML uses the biased gamma and the standard statistic, so the `UG`/`RLS`
suffixes do not apply. The default `fiml.convention = "observed"` uses
observed saturated and model information throughout. Use
`fiml.convention = "lavaan"` to reproduce lavaan 0.7-2’s inspected
robust-test spectrum; the convention is printed with the result. The
same option is available in `pvalues_nested()`, where
`A.method = "delta"` is the default restriction map.

``` r
HS <- lavaan::HolzingerSwineford1939
set.seed(313)
HS$x1[sample(nrow(HS), 60)] <- NA
fit_fiml <- lavaan::cfa(model, HS, missing = "fiml", estimator = "MLR")
pvalues(fit_fiml)
#>     peba4_ml 
#> 3.301168e-07
#> estimator: ML (MLR) (FIML) | data: continuous | information: observed | df: 24 | FIML convention: observed
```

For a single FIML model, the lavaan convention reproduces
`lavInspect(fit_fiml, "UGamma")`. It does not necessarily reproduce the
scalar Yuan–Bentler–Mplus test stored by a default MLR fit. For nested
FIML models, the lavaan convention reproduces the public Satorra-2000
difference test.

Categorical and mixed-indicator models use lavaan’s inspected `UGamma`
directly. The DWLS, ULS, and full-WLS estimator families are supported
experimentally, including single- and multigroup models and listwise or
pairwise missingness. Pairwise support reproduces lavaan’s pairwise
sample-statistic calculation; it is not FIML and is not a claim of
generally MAR-valid inference.

``` r
HSord <- lavaan::HolzingerSwineford1939
ordered_names <- paste0("x", 1:9)
HSord[ordered_names] <- lapply(
  HSord[ordered_names],
  function(x) ordered(cut(x, 3))
)
fit_ordinal <- lavaan::cfa(model, HSord, ordered = ordered_names)
pvalues(fit_ordinal)
#>     peba4_ml
#> 7.156389e-06
#> estimator: DWLS (WLSMV) | data: categorical | information: expected | df: 24
```

Nested categorical models support the Satorra-2000 delta-map
construction:
`pvalues_nested(m0, m1, method = "2000", A.method = "delta")`.

## References

Foldnes, N., Moss, J., & Grønneberg, S. (2024). Improved goodness of fit
procedures for structural equation models. Structural Equation Modeling:
A Multidisciplinary Journal, 1–13.
<https://doi.org/10.1080/10705511.2024.2372028>

Foldnes, N., Grønneberg, S., & Moss, J. (2026). Penalized eigenvalue
block averaging: Extension to nested model comparison and Monte Carlo
evaluations. Behavior Research Methods.
<https://doi.org/10.3758/s13428-026-02968-4>

Foldnes, N., & Grønneberg, S. (2018). Approximating Test Statistics
Using Eigenvalue Block Averaging. Structural Equation Modeling, 25(1),
101–114. <https://doi.org/10.1080/10705511.2017.1373021>

Du, H., & Bentler, P. M. (2022). 40-Year Old Unbiased Distribution Free
Estimator Reliably Improves SEM Statistics for Nonnormal Data.
Structural Equation Modeling: A Multidisciplinary Journal, 29(6),
872–887. <https://doi.org/10.1080/10705511.2022.2063870>

Grønneberg, S., & Foldnes, N. (2019). Testing Model Fit by Bootstrap
Selection. Structural Equation Modeling, 26(2), 182–190.
<https://doi.org/10.1080/10705511.2018.1503543>

Marcoulides, K. M., Foldnes, N., & Grønneberg, S. (2020). Assessing
Model Fit in Structural Equation Modeling Using Appropriate Test
Statistics. Structural Equation Modeling, 27(3), 369–379.
<https://doi.org/10.1080/10705511.2019.1647785>

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. Journal of Statistical Software, 48(2), 1–36.
<https://doi.org/10.18637/jss.v048.i02>

## How to Contribute or Get Help

If you encounter a bug, have a feature request or need some help, open a
[Github issue](https://github.com/JonasMoss/semTests/issues). Create a
pull requests to contribute.
