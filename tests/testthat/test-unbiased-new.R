values_nested_groups <- list("2000_SB_RLS" = pvalues_nested(m0, m1, "2000", "SB_RLS"),
     "2000_SB_UG_RLS" = pvalues_nested(m0, m1, "2000", "SB_UG_RLS"))

values_nested_no_groups <- list("2000_SB_RLS" = pvalues_nested(m0_no_groups, m1_no_groups, "2000", "SB_RLS"),
                             "2000_SB_UG_RLS" = pvalues_nested(m0_no_groups, m1_no_groups, "2000", "SB_UG_RLS"))

values_no_groups <- list("SB_RLS" = pvalues(object, "SB_RLS"),
               "SB_UG_RLS" = pvalues(object, "SB_UG_RLS"))

values_groups <- list("SB_RLS" = pvalues(m0, "SB_RLS"),
                      "SB_UG_RLS" = pvalues(m0, "SB_UG_RLS"))

#saveRDS(values_nested_no_groups, "values_nested_no_groups.Rds")
#saveRDS(values_nested_groups, "values_nested_groups.Rds")
#saveRDS(values_no_groups, "values_no_groups.Rds")
#saveRDS(values_groups, "values_groups.Rds")

without_provenance <- function(x) {
  lapply(x, function(value) {
    attr(value, "semtests") <- NULL
    value
  })
}

test_that("testhat", {
  expect_equal(
    without_provenance(readRDS("values_no_groups.Rds")),
    without_provenance(values_no_groups)
  )
  expect_equal(
    without_provenance(readRDS("values_groups.Rds")),
    without_provenance(values_groups)
  )
  expect_equal(
    without_provenance(readRDS("values_nested_no_groups.Rds")[names(values_nested_no_groups)]),
    without_provenance(values_nested_no_groups)
  )
  expect_equal(
    without_provenance(readRDS("values_nested_groups.Rds")[names(values_nested_groups)]),
    without_provenance(values_nested_groups)
  )
})
