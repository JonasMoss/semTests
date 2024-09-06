hs_model_0 <- " visual  =~ x1 + a*x2 + x3
                textual =~ x4 + a*x5 + a*x6
                speed   =~ x7 + x8 + x9 "

hs_model_1 <- " visual  =~ x1 + x2 + x3
                textual =~  x4 + x5 + x6
                speed   =~  x7 + x8 + x9 "


## Estimation that IS allowed.
m1 <- lavaan::cfa(hs_model_1,
                  data = lavaan::HolzingerSwineford1939, estimator = "MLM"
)

m0 <- lavaan::cfa(hs_model_0,
                  data = lavaan::HolzingerSwineford1939, estimator = "MLM"
)


lavaan::anova(m1, m0, method = "satorra.2000", scaled.shifted = FALSE)
pvalues_nested(m0, m1, tests = "SB")

pvalues_nested(m0, m1, tests = "SB", method = "2001")
lavaan::anova(m1, m0, method = "satorra.bentler.2001")
