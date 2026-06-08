## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on macos (release), windows (release), and ubuntu (devel, release, oldrel).

## Submission notes

This is a maintenance release. It fixes compatibility with the upcoming
lavaan 0.7-1, in which `lavaan::lavTest()` returns a named list of test
results; the chi-square statistic is now extracted in a way that works with
both old and new lavaan versions.
