## R CMD check results

0 errors | 0 warnings | 0 notes

Tested with GitHub Actions on macos (release), windows (release), and ubuntu
(devel, release, oldrel), and locally on ubuntu.

## Submission notes

* Ensures the forthcoming lavaan release can be published (compatibility fix).
* Adds new features (more estimators, plus FIML for missing data).
* Removes dependencies (Imports is now just lavaan and methods).

See NEWS.md for details.
