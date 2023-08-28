library("semselector")

## One model, no restrictions
set.seed(313)
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

unbiased <- ugamma_non_nested(object, unbiased = TRUE)
biased <- ugamma_non_nested(object, unbiased = FALSE)
mean(abs((unbiased - biased)/biased))

## Two models

f0 <- ' visual  =~ a*x1 + a*x2 + x3
        textual =~ b*x4 + b*x5 + x6
        speed   =~ x7 + x8 + x9 '

f1 <- ' visual  =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6
        speed   =~ x7 + x8 + x9 '

m0 <- lavaan::cfa(f0, data = lavaan::HolzingerSwineford1939)
m1 <- lavaan::cfa(f1, data = lavaan::HolzingerSwineford1939)

ugamma_nested(m0, m1)

## Try two
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

ubiased <- ugamma_nested(m0, m1)
uunbiased <- ugamma_nested(m0, m1, unbiased = TRUE)

