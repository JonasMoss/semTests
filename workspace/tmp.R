hs_model <- " visual  =~ x1 + a*x2 + x3
              textual =~ x4 + x5 + a*x6
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

hs_model0 <- " visual  =~ x1 + a*x2 + x3
              textual =~ x4 + b*x5 + a*x6
              speed   =~ x7 + b*x8 + x9 "

hs_model1 <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "


m1 <- lavaan::cfa(hs_model1,
                  data = lavaan::HolzingerSwineford1939,
                  estimator = "MLM"
)

m0 <- lavaan::cfa(hs_model0,
                  data = lavaan::HolzingerSwineford1939,
                  estimator = "MLM"
)


sum(diag(lav_ugamma_nested(m0, m1)))
lavaan:::lav_test_diff_Satorra2000(m1, m0)$trace.UGamma


a = lavaan:::lav_test_diff_satorra(m0, m1)

(124.04 - 115.85) / 7.7257
