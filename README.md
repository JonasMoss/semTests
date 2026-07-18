
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semTests <img src="man/figures/logo.png" align="right" width="210" height="130" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/semTests)](https://cran.r-project.org/package=semTests)
[![R-CMD-check](https://github.com/JonasMoss/semTests/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JonasMoss/semTests/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/JonasMoss/semTests/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JonasMoss/semTests?branch=main)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

An R package for robust tests of overall fit and nested model
comparisons in structural equation models. It works with `lavaan` and
implements the penalized eigenvalue block averaging and penalized
regression *p* values of Foldnes, Moss, and Grønneberg (2025), along
with their extension to nested model comparison in Foldnes, Grønneberg,
and Moss (2026).

**News.** You can now use `semTests` with GLS and ULS, categorical DWLS,
ULS, or WLS models, and continuous FIML fits. Nested comparisons are
included too. See `?semTests-support` for the details.

## Installation

Install the released version from CRAN:

``` r
install.packages("semTests")
```

Or try the development version:

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/semTests")
```

## Examples

Here is a quick CFA:

``` r
library("semTests")
model <- "visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9"
object <- lavaan::cfa(model, lavaan::HolzingerSwineford1939, estimator = "MLM")
pvalues(object, tests = c("SB", "SS", "PEBA4"))
#>       sb_rls       ss_rls    peba4_rls
#> 1.726558e-07 9.251198e-07 5.165737e-07
#> estimator: ML (MLM) | data: continuous | information: expected | df: 24
```

For a nested comparison, pass the constrained model first:

``` r
constrained <- "visual  =~ x1 + x2 + x3
                textual =~ a*x4 + a*x5 + a*x6
                speed   =~ x7 + x8 + x9"
m1 <- lavaan::cfa(model, lavaan::HolzingerSwineford1939, estimator = "MLM")
m0 <- lavaan::cfa(constrained, lavaan::HolzingerSwineford1939, estimator = "MLM")
pvalues_nested(m0, m1, tests = c("SB", "SS", "PEBA2"))
#>     sb_rls     ss_rls  peba2_rls
#> 0.01023556 0.01092401 0.01051737
#> estimator: ML (MLM) | data: continuous | information: expected | df: 2 | nested (method 2000)
```

The fits need to use the same data, estimator, and fitting choices. If
they arrive in the wrong order, `semTests` lets you know and swaps them.

Missing values are fine too. Fit the model with `missing = "fiml"` and
call `pvalues()` as usual. `semTests` uses observed information by
default, following Savalei (2010). The [FIML
vignette](vignettes/fiml-missing-data.Rmd) explains this choice and
shows the other options.

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

Ordered indicators work in much the same way:

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

## More examples

- `vignette("continuous-data", package = "semTests")`
- `vignette("categorical-data", package = "semTests")`
- `vignette("fiml-missing-data", package = "semTests")`
- `vignette("semTests", package = "semTests")`

The help page `?semTests-support` gives the exact supported
combinations.

## References

Foldnes, N., Moss, J., & Grønneberg, S. (2025). Improved goodness of fit
procedures for structural equation models. Structural Equation Modeling:
A Multidisciplinary Journal, 32(1), 1–13.
<https://doi.org/10.1080/10705511.2024.2372028>

Foldnes, N., Grønneberg, S., & Moss, J. (2026). Penalized eigenvalue
block averaging: Extension to nested model comparison and Monte Carlo
evaluations. Behavior Research Methods, 58, article 107.
<https://doi.org/10.3758/s13428-026-02968-4>

Savalei, V. (2010). Expected versus observed information in SEM with
incomplete normal and nonnormal data. Psychological Methods, 15(4),
352–367. <https://doi.org/10.1037/a0020143>

## How to Contribute or Get Help

Found a bug or have an idea? Open a [GitHub
issue](https://github.com/JonasMoss/semTests/issues). Pull requests are
welcome too.
