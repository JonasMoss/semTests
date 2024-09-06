model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10],
                      scaled.test = "browne.residual.nt.model")


pvalues(object, "SB_RLS")
pvalues(object, "SB")

pvalues(object, "SS_RLS")

pvalues(object, "SS")


pvalues(object, c("SB_RLS", "SB"))

























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


##
## Verify gamma.
##

unbiased = FALSE
lavoptions = lavaan::lavInspect(m1, "options")
lavdata = m1@Data
lavoptions$gamma.unbiased <- unbiased
gamma_list <- lavaan:::lav_samplestats_from_data(lavdata, lavoptions)@NACOV
x3 <- as.matrix(Matrix::bdiag(gamma_list))
sum(abs(x2 - x3))

x1 <- get_gamma(m1, FALSE)
x2 <- lavaan:::lav_matrix_bdiag(lavaan::lavInspect(m1, "gamma"))
sum(abs(x1-x2))


##
## Verify U.
##

# Compute scaling factor
fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

gamma <- lavaan:::lav_samplestats_from_data(lavdata, lavoptions)@NACOV
gamma_rescaled <- gamma
for (i in (seq_along(gamma))) {
  gamma_rescaled[[i]] <- fg[i] * gamma_rescaled[[i]]
}

gamma <- as.matrix(Matrix::bdiag(gamma_rescaled))


ug1 <- lavaan::lavInspect(m1, "U") %*% gamma
ug2 <- lavaan::lavInspect(m1, "UGamma")
test <- \(a, b) sum(abs(a - b))
test(ug1, ug2)


df <- length(lambdas)
tr_ug <- sum(lambdas)
tr_ug2 <- sum(lambdas^2)
a <- sqrt(df / tr_ug2)
b <- df - sqrt(df * tr_ug^2 / tr_ug2)
t3 <- unname(chisq * a + b)
1 - stats::pchisq(t3, df = df)
