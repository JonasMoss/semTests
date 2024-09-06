tester <- \(object) {
  gammas <- get_gamma(object, collapse = FALSE)
  gamma1 <- gamma_to_gamma_unbiased(gammas, object)
  gamma2 <- get_gamma(object, unbiased = TRUE, collapse = FALSE)

  sum(abs(gamma1[[1]] - gamma2[[1]])) + sum(abs(gamma1[[2]] - gamma2[[2]]))
}

hs_model <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "

m1 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939,
                  group = "school", estimator = "MLM"
)

m0 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939,
                  group = "school", estimator = "ML", group.equal = "loadings", test = "satorra.bentler"
)

m3 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939, estimator = "ML", test = "satorra.bentler"
)


tester(m0)
tester(m1)
tester(m3)



gammas <- get_gamma(object, collapse = FALSE)
gamma1 <- gamma_to_gamma_unbiased(gammas, object)
gamma2 <- get_gamma(object, unbiased = TRUE, collapse = FALSE)

sum(abs(gamma1[[1]] - gamma2[[1]]))
sum(abs(gamma1[[2]] - gamma2[[2]]))
