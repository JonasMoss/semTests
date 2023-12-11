
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semTests <img src="man/figures/logo.png" align="right" width="210" height="130" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/semTests)](https://cran.r-project.org/package=semTests)
[![R-CMD-check](https://github.com/JonasMoss/semTests/workflows/R-CMD-check/badge.svg)](https://github.com/JonasMoss/semTests/actions)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

An R package for goodness of fit testing of structural equation models.
Built on top of `lavaan`.

## Installation

Use the following command from inside `R`:

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/semTests")
```

## Usage

Call the `library` function, create a `lavaan` model, and run the
`pvalues` function.

``` r
library("semTests")
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 200
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")
pvalues(object)
#> peba_half2_trad peba_half4_trad  peba_half2_rls  peba_half4_rls 
#>      0.04182486      0.04370615      0.04235143      0.04424139
```

## References

Foldnes, N., & Grønneberg, S. (2018). Approximating Test Statistics
Using Eigenvalue Block Averaging. Structural Equation Modeling, 25(1),
101–114. <https://doi.org/10.1080/10705511.2017.1373021>

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
pull requests to contribute. This project follows a [Contributor Code of
Conduct](https://www.contributor-covenant.org/version/1/4/code-of-conduct.md).
