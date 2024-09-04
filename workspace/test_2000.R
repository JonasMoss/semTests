hs_model_1 <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "

hs_model_0 <- " visual  =~ x1 + a*x2 + x3
              textual =~ x4 + a*x5 + a*x6
              speed   =~ x7 + x8 + x9 "

m1 <- lavaan::cfa(hs_model_1,
                  data = lavaan::HolzingerSwineford1939, test = "SB"
)

m0 <- lavaan::cfa(hs_model_0,
                  data = lavaan::HolzingerSwineford1939, test = "SB"
)

hs_model <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "

m1 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939,
                  group = "school", estimator = "MLM"
)

m0 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939,
                  group = "school", estimator = "MLM", group.equal = "loadings"
)


pvalues_nested(m0, m1, tests = "SB", method = "2000")
lavaan::anova(m1, m0, method = "satorra.2000", scaled.shifted = FALSE)

ug <- ugamma_nested(m0, m1, "2000", 1)[[1]]
trace <- \(x) sum(diag(x))

scale <- attr(lavaan::anova(m1, m0, method = "satorra.2000", scaled.shifted = FALSE), "scale")[2]
x <- (lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df"))/scale

x <- trace(ug)^2 / trace(ug%*%ug) /scale

trace(ug) / x

##
## Verify scaled and shifted.
##

chisq <- make_chisqs("ml", m0, m1)
values <- Re(eigen(lav_ugamma_nested_2000(m0, m1, get_gamma(m0, FALSE, collapse = FALSE)))$values)
lambdas <- values[1:(lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df"))]
scaled_and_shifted(chisq,lambdas)

lavaan::anova(m1, m0, method = "satorra.2000", scaled.shifted = TRUE)


df <- length(lambdas)
tr_ug <- sum(lambdas)
tr_ug2 <- sum(lambdas^2)
a <- sqrt(df / tr_ug2)
b <- df - sqrt(df * tr_ug^2 / tr_ug2)
t3 <- unname(chisq * a + b)
1 - stats::pchisq(t3, df = df)
