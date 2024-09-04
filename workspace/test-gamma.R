### =======================================
### Computed gammas are not already scaled.
### =======================================


x <- lavaan::HolzingerSwineford1939
g1 <- x[x$school == "Pasteur", 7:15]
g2 <- x[x$school == "Grant-White", 7:15]

gamma1 <- gamma_est_adf(as.matrix(g1))
gamma2 <- gamma_est_adf(as.matrix(g2))

hs_model <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "

m1 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939,
                  group = "school", estimator = "ML", test = "satorra.bentler"
)

gamma_all <- lavInspect(m1, "gamma")

gamma_all[[2]][10:54, 10:54] / gamma2

### =======================================
### Computed gammas are not already scaled.
### =======================================


gamma_1 <- lav_ugamma_nested_2000(m0, m1, get_gamma(m0, FALSE, FALSE), strange = FALSE)
gamma_2 <- lav_ugamma_nested_2000(m0, m1, get_gamma(m0, FALSE, FALSE))

tr(gamma_1)
tr(gamma_2)
sum(abs(gamma_1 - gamma_2))
