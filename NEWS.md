# semTests 1.0.0

A major release: semTests now has a self-contained numerical core and a minimal
dependency footprint -- `Imports` is just `lavaan` and the base `methods`
package.

* Fixed compatibility with lavaan 0.7-1, in which `lavaan::lavTest()` returns a
  named list of test results rather than a single flat object (the reweighted
  least squares chi-square statistic was previously read as `NULL`). The
  statistic is now extracted in a way that works with both old and new lavaan
  return structures.
* The eigenvalue-based p-values no longer depend on CompQuadForm. The upper tail
  of the quadratic form `sum lambda_j chi^2_1` is computed internally by
  `imhof_pvalue()`, which is more accurate in the tail, where CompQuadForm
  degrades (its `imhof` returns quadrature noise, including negative
  probabilities, once the p-value drops below roughly `1e-6`, and `davies`
  underflows to zero). Positive spectra use Ruben's series (exact to machine
  precision, including the tail); mixed-sign spectra use the Imhof integral in
  the body with a Lugannani-Rice saddlepoint fallback in the tail.
* Nested tests now read the reference spectrum from the reduced `m x m`
  restriction-space matrix instead of the full `q x q` UGamma, where `m` is the
  number of restrictions. The result is identical to machine precision but
  faster on large models, and it removes the need for `RSpectra`.
* Replaced the internal use of `MASS::ginv` with a self-contained
  `generalized_inverse()`.
* Removed the unused internal bootstrap helpers.
* Dropped the `CompQuadForm`, `Matrix`, `MASS`, `RSpectra`, `progressr`, and
  `future.apply` dependencies, plus the `psych` suggestion (examples now use the
  built-in `HolzingerSwineford1939` data). `Imports` is now `lavaan` and
  `methods`.
* Documentation fixes and a substantially expanded test suite.
