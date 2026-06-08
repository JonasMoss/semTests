# semTests 0.7.2

* Fixed compatibility with the upcoming lavaan 0.7-1, in which
  `lavaan::lavTest()` returns a named list of test results rather than a single
  flat object. This previously caused the reweighted least squares chi-square
  statistic to be read as `NULL`. The statistic is now extracted in a way that
  works with both old and new lavaan return structures.
* Removed the unused internal bootstrap helpers, which dropped the `progressr`
  and `future.apply` dependencies.
* Minor documentation and code cleanup.
