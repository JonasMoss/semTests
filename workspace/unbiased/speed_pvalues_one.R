pvalues_one_with <- function(object) {
  chisq <- lavaan::fitmeasures(object, "chisq")
  ug <- ugamma_non_nested(object)
  ug_unbiased <- ugamma_non_nested(object, TRUE)
  df <- lavaan::fitmeasures(object, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]
  eigenps <- eigen_pvalues(chisq, lambdas)
  lambdas_unbiased <- Re(eigen(ug_unbiased)$values)[seq(df)]
  eigenps_unbiased <- eigen_pvalues(chisq, lambdas_unbiased)

  c(
    pstd = unname(lavaan::fitmeasures(object, "pvalue")),
    psb = eigenps$psb,
    phalf = eigenps$phalf,
    pfull = eigenps$pfull,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(chisq, lambdas),
    psb_unb = eigenps_unbiased$psb,
    phalf_unb = eigenps_unbiased$phalf,
    pfull_unb = eigenps_unbiased$pfull,
    psf_unb = scaled_f(chisq, lambdas_unbiased),
    pss_unb = scaled_and_shifted(chisq, lambdas_unbiased)
  )
}

pvalues_one_without <- function(object) {
  chisq <- lavaan::fitmeasures(object, "chisq")
  ug <- ugamma_non_nested(object)
  df <- lavaan::fitmeasures(object, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]
  eigenps <- eigen_pvalues(chisq, lambdas)
  c(
    pstd = unname(lavaan::fitmeasures(object, "pvalue")),
    psb = eigenps$psb,
    phalf = eigenps$phalf,
    pfull = eigenps$pfull,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(chisq, lambdas)
  )
}

n <- 50
nn <- sample(nrow(psych::bfi), n)
object <- lavaan::sem(model, psych::bfi[nn, 1:10], estimator = "MLM")

microbenchmark::microbenchmark(pvalues_one_with(object), pvalues_one_without(object))



object <- lavaan::sem(model, lavaan::simulateData(model, sample.nobs = 50), estimator = "MLM")

ug <- ugamma_non_nested(object)
ug_unbiased <- ugamma_non_nested(object, TRUE)
lambdas <- Re(eigen(ug)$values)[seq(df)]
lambdas_unbiased <- Re(eigen(ug_unbiased)$values)[seq(df)]

plot(lambdas_unbiased, col = "blue")
points(lambdas)
lines(big_lambdas)
abline(h = mean(big_lambdas))



object_big <- lavaan::sem(model,
                      lavaan::simulateData(model, sample.nobs = 100000, kurtosis = 22),
                      estimator = "MLM")
ug_big <- ugamma_non_nested(object_big)
lambdas_big <- Re(eigen(ug)$values)[seq(df)]


object <- lavaan::sem(model,
                          lavaan::simulateData(model, sample.nobs = 50, kurtosis = 22),
                          estimator = "MLM")
ug <- ugamma_non_nested(object)
ug_unbiased <- ugamma_non_nested(object, TRUE)
lambdas <- Re(eigen(ug)$values)[seq(df)]
lambdas_unbiased <- Re(eigen(ug_unbiased)$values)[seq(df)]

plot(lambdas_unbiased, col = "blue")
points(lambdas)
lines(lambdas_big)


library("future.apply")
library("quadagree")
plan(multisession, workers = availableCores() - 2)

results <- future.apply::future_replicate(1000, {
  object <- lavaan::sem(model,
                        lavaan::simulateData(model, sample.nobs = 50, kurtosis = 22),
                        estimator = "MLM")
  lavaan::fitmeasures(object, "chisq")
})

df <- 34
hist(results, freq = FALSE, breaks = 100)
y <- seq(5, 100, by = 0.5)
lines(y, dchisq(y, 34))
