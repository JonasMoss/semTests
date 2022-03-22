
<!-- README.md is generated from README.Rmd. Please edit that file -->

    #> Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    #> variances are negative

    #> Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    #> variances are negative

# semselector <img src="man/figures/logo.png" align="right" width="210" height="130" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/semselector)](https://cran.r-project.org/package=semselector)
[![R-CMD-check](https://github.com/JonasMoss/semselector/workflows/R-CMD-check/badge.svg)](https://github.com/JonasMoss/semselector/actions)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

An R package for goodness of fit testing for structural equation models

## Installation

Use the following command from inside `R`:

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/semselector")
```

## Usage

Call the `library` function, prepare data into matrix form, and run the
`semselector` function.

``` r
library("semselector")
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 200
object <- lavaan::sem(model, psych::bfi[1:n, 1:10])
pvalues(object)
#>        pml        psb      pfull      phalf        pcf        pss 
#> 0.01038449 0.03688981 0.06563318 0.05535394 0.06540780 0.06287733
```

You can calculate the distribution of *p*-values using the following
command:

``` r
library("semselector")
set.seed(313)
boots <- bootstrapper(object, n_reps = 5000)

library("ggplot2")
theme_set(theme_minimal())
ggplot(reshape2::melt(as.data.frame(t(boots))), aes(value)) +
  geom_histogram() +
  xlab("p-value") +
  facet_wrap(~variable, nrow = 3)
```

    #> No id variables; using all as measure variables
    #> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

<img src="man/figures/README-plot_eval-1.png" width="750px" />

## References

-   **fill in**.

## How to Contribute or Get Help

If you encounter a bug, have a feature request or need some help, open a
[Github issue](https://github.com/JonasMoss/semselector/issues). Create
a pull requests to contribute. This project follows a [Contributor Code
of
Conduct](https://www.contributor-covenant.org/version/1/4/code-of-conduct.md).
