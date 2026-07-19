hs_model <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "

## Estimation that IS allowed.
m1 <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "MLM"
)

m0 <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "MLM", group.equal = "loadings"
)

object <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  estimator = "MLM"
)


## Estimation that isn't allowed.
m1_ <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "GLS"
)

m0_ <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "GLS", group.equal = "loadings"
)



## Nested without groups

hs_model_no_groups <- " visual  =~ x1 + a*x2 + x3
              textual =~ x4 + a*x5 + a*x6
              speed   =~ x7 + a*x8 + x9 "

m1_no_groups <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939, estimator = "MLM"
)

m0_no_groups <- lavaan::cfa(hs_model_no_groups,
  data = lavaan::HolzingerSwineford1939, estimator = "MLM"
)
