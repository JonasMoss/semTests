## R CMD check results

0 errors | 0 warnings | 1 note

Tested on macos (release), windows (release), and ubuntu (devel, release, oldrel).

* This is a new release.
* The note "Unexported objects imported by ':::' calls:
     'lavaan:::lav_model_information' 'lavaan:::lav_object_gamma'
     'lavaan:::lav_test_diff_A'
     See the note in ?`:::` about the use of this operator."
  is not avoidable; we have attempted to use only exported functions,
  but the functionality we need is not exported.
