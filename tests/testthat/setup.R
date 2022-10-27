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

object <- m0
data <- bollen_stine_transform(m0)


## Estimation that isn't allowed.
m1_ <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "GLS"
)

m0_ <- lavaan::cfa(hs_model,
  data = lavaan::HolzingerSwineford1939,
  group = "school", estimator = "GLS", group.equal = "loadings"
)
