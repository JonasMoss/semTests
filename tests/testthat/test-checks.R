m <- "visual =~ x1 + x2 + x3
      textual =~ x4 + x5 + x6
      speed =~ x7 + x8 + x9"

test_that("warn_fiml_information warns only for FIML with expected information", {
  HS <- lavaan::HolzingerSwineford1939
  set.seed(1)
  HS$x1[sample(nrow(HS), 60)] <- NA

  fiml_observed <- lavaan::cfa(m, HS, missing = "fiml", estimator = "MLR")
  fiml_expected <- lavaan::cfa(m, HS, missing = "fiml", estimator = "ML",
                               information = "expected")
  complete      <- lavaan::cfa(m, na.omit(HS), estimator = "ML")  # expected, but complete

  expect_warning(warn_fiml_information(fiml_expected), "MCAR")
  expect_silent(warn_fiml_information(fiml_observed))   # FIML but observed -> fine
  expect_silent(warn_fiml_information(complete))        # expected but complete -> fine
})
