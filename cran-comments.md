## Submission

semTests 2.0.0 adds score and criterion-difference tests (`semtest()`) and
gives nested tests a default restriction matrix that does not depend on the
units of the variables. It requires lavaan (>= 0.7-3), deprecates the
`fiml.convention` argument, and fixes three bugs in nested comparisons.

## Test environments

* Local: Ubuntu 26.04, R 4.6.1
* GitHub Actions: macOS (arm64), Windows, Ubuntu (R-devel, release, oldrel),
  and a noLD build

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

None.
