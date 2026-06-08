# semTests 0.7.3

* The eigenvalue-based p-values no longer depend on CompQuadForm. The upper
  tail of the quadratic form `sum lambda_j chi^2_1` is now computed internally
  by `imhof_pvalue()`. It is exact in the tail, where CompQuadForm degrades:
  its `imhof` returns quadrature noise (including negative probabilities) once
  the p-value drops below roughly `1e-6`, and `davies` underflows to zero. For
  positive eigenvalues `imhof_pvalue()` uses Ruben's series (exact to machine
  precision, including the tail); mixed-sign spectra use the Imhof integral in
  the body with a Lugannani-Rice saddlepoint fallback in the tail.
* Dropped the `CompQuadForm` and `Matrix` dependencies.

# semTests 0.7.2

* Fixed compatibility with the upcoming lavaan 0.7-1, in which
  `lavaan::lavTest()` returns a named list of test results rather than a single
  flat object. This previously caused the reweighted least squares chi-square
  statistic to be read as `NULL`. The statistic is now extracted in a way that
  works with both old and new lavaan return structures.
* Removed the unused internal bootstrap helpers, which dropped the `progressr`
  and `future.apply` dependencies.
* Minor documentation and code cleanup.
