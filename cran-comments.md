## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on macos (release), windows (release), and ubuntu (devel, release, oldrel).

## Submission notes

This release drops the CompQuadForm and Matrix dependencies. The upper-tail
probabilities of the quadratic forms underlying the eigenvalue-based p-values
are now computed internally, which is also more accurate in the tail, where
CompQuadForm's routines lose precision.
