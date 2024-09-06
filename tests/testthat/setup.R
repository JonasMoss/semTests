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

model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

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
