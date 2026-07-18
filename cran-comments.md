## Resubmission

This is a resubmission of a package that was archived on CRAN on 2026-07-15
("issues were not corrected in time"). Version 1.0.0 restores compatibility with
the current CRAN lavaan (>= 0.7-2) and removes a stale, now-invalid DOI from the
Browne (1974) reference.

## R CMD check results

Local `R CMD check --as-cran` (Ubuntu, R 4.5.2) reports no ERRORs, WARNINGs, or
NOTEs beyond the expected incoming-feasibility note that the package is a new
submission and was previously archived on CRAN.

## Notes

* Requires lavaan (>= 0.7-2), now available on CRAN.
* The GLS, ULS, categorical DWLS/ULS, and single-group FIML estimators (single
  and nested) are now documented as stable, having been validated to numerical
  tolerance against an independent implementation.
* Full WLS/ADF and nested method "2001" are refused pending further work, each
  with a message pointing to the supported alternative.

See NEWS.md for the full list of changes.
