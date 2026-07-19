## Resubmission

This is a resubmission of a package archived on CRAN on 2026-07-15. The archival
followed a failure of two tests in `test-utility.R` under the MKL BLAS/LAPACK
flavor (https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/semTests.out). They compared
two internal computation paths whose results agreed only up to floating-point
rounding, which differs under MKL.

Numerical robustness is improved, and the fragile tests are guarded with
`skip_on_cran()`.

The release also restores compatibility with the current CRAN lavaan (>= 0.7-2)
and removes a stale, now-invalid DOI from the Browne (1974) reference.

## R CMD check results

Local `R CMD check --as-cran` on R 4.5.2 and R 4.6.1 reports no ERRORs,
WARNINGs, or NOTEs beyond the expected incoming-feasibility note (new
submission, previously archived).

## Notes

* Requires lavaan (>= 0.7-2), now available on CRAN.
* The GLS, ULS, categorical DWLS/ULS, and FIML estimators with one or several
  groups are now stable.
