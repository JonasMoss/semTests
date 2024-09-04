hs_model <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "


m1_ml <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "ML"
)

m1_ml_sb <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "ML", test = "satorra.bentler"
)

m1_mlm <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "MLM"
)


m0_ml <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "ML", group.equal = "loadings"
)


m0_mlm <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "MLM", group.equal = "loadings"
)


m0_ml_sb <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "ML", group.equal = "loadings", test = "satorra.bentler"
)


object_ml <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939, estimator = "ML"
)

object_ml_sb <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939, estimator = "ML", test = "satorra.bentler"
)

object_mlm <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939, estimator = "MLM"
)




expect_error(get_gamma(m1_ml, FALSE, FALSE, m0_ml))

expect_equal(
  sum(get_gamma(m1_ml, FALSE, FALSE, m0_mlm)[[1]]),
  sum(get_gamma(m1_ml, FALSE, FALSE, m0_ml_sb)[[1]]))

expect_equal(
  sum(get_gamma(m1_mlm, FALSE, FALSE, m0_ml)[[1]]),
  sum(get_gamma(m1_mlm, FALSE, FALSE, m0_mlm)[[1]])
)

expect_equal(
  sum(get_gamma(m1_mlm, FALSE, FALSE, m0_ml_sb)[[1]]),
  sum(get_gamma(m1_ml_sb, FALSE, FALSE, m0_ml)[[1]])
)

expect_equal(
  sum(get_gamma(m1_ml_sb, FALSE, FALSE, m0_mlm)[[1]]),
  sum(get_gamma(m1_ml_sb, FALSE, FALSE, m0_ml_sb)[[1]])
)

expect_error(get_gamma(m1_ml, FALSE, TRUE))
expect_equal(sum(get_gamma(m1_mlm, FALSE, TRUE)), sum(get_gamma(m1_ml_sb, FALSE, TRUE)))

expect_error(get_gamma(object_ml, FALSE, TRUE))
expect_equal(sum(get_gamma(object_mlm, FALSE, TRUE)), sum(get_gamma(object_ml_sb, FALSE, TRUE)))
