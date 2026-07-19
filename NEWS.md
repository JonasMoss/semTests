# semTests 1.0.0

Earlier versions marked several estimators as experimental. This release
promotes them to stable.

## New features

* GLS, ULS, categorical DWLS and ULS, mixed continuous and ordinal indicators,
  and FIML, each with one or several groups, now work for both single models and
  nested comparisons.
* A `fiml.convention` argument selects the FIML information convention. The
  default disagrees with lavaan but has better statistical properties; set
  `fiml.convention = "lavaan"` to reproduce lavaan's numbers exactly.
* New vignettes cover measurement invariance with ordinal indicators and latent
  growth models under FIML.

## Bug fixes

* The reweighted least-squares statistic was read as `NULL` under lavaan 0.7-2
  and is now correct. semTests requires `lavaan (>= 0.7-2)`.

## Removed

* Nested method 2001 is withdrawn because of its poor performance.
* Full WLS is refused because its correction reduces to the ordinary chi-square.
