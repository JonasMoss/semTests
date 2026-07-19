# semTests 1.0.0

## New features

* GLS, ULS, categorical DWLS and ULS, mixed indicators, and FIML with one or
  several groups now work for single models and nested comparisons.
* A `fiml.convention` argument selects the FIML information convention, with
  default disagreeing with lavaan, but with better statistical properties.
* New vignettes cover measurement invariance with ordinal indicators and latent
  growth models under FIML.

## Bug fixes

* The reweighted least-squares statistic was read as `NULL` under lavaan 0.7-2
  and is now correct. semTests requires `lavaan (>= 0.7-2)`.

## Removed

* Nested method 2001 is withdrawn because of its poor performance.
* Full WLS is refused because its correction reduces to the ordinary chi-square.
