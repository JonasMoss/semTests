pvalues_two_with <- function(m0, m1) {
  if (m0@Options$estimator != "ML" || m1@Options$estimator != "ML" ||
      m0@Options$se == "standard" || m1@Options$se == "standard") {
    stop("Only the 'ML' estimator has currently tested.")
  }
  chisq <- lavaan::fitmeasures(m0, "chisq") - lavaan::fitmeasures(m1, "chisq")
  ug <- ugamma_nested(m0, m1)
  ug_unbiased <- ugamma_nested(m0, m1, TRUE)
  df <- lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]
  eigenps <- eigen_pvalues(chisq, lambdas)
  lambdas_unbiased <- Re(eigen(ug_unbiased)$values)[seq(df)]
  eigenps_unbiased <- eigen_pvalues(chisq, lambdas)

  c(
    pstd = unname(1 - stats::pchisq(chisq, df)),
    psb = eigenps$psb,
    phalf = eigenps$phalf,
    pfull = eigenps$pfull,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(chisq, lambdas),
    psb_unb  = eigenps_unbiased$psb,
    phalf_unb = eigenps_unbiased$phalf,
    pfull_unb = eigenps_unbiased$pfull,
    psf_unb  = scaled_f(chisq, lambdas_unbiased),
    pss_unb  = scaled_and_shifted(chisq, lambdas_unbiased)
  )
}

pvalues_two_without <- function(m0, m1) {
  if (m0@Options$estimator != "ML" || m1@Options$estimator != "ML" ||
      m0@Options$se == "standard" || m1@Options$se == "standard") {
    stop("Only the 'ML' estimator has currently tested.")
  }
  chisq <- lavaan::fitmeasures(m0, "chisq") - lavaan::fitmeasures(m1, "chisq")
  ug <- ugamma_nested(m0, m1)
  ug_unbiased <- ugamma_nested(m0, m1, TRUE)
  df <- lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]
  eigenps <- eigen_pvalues(chisq, lambdas)

  c(
    pstd = unname(1 - stats::pchisq(chisq, df)),
    psb = eigenps$psb,
    phalf = eigenps$phalf,
    pfull = eigenps$pfull,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(chisq, lambdas)
  )
}

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

pvalues_two_with(m0, m1)

microbenchmark::microbenchmark(pvalues_one_with(object),
                               pvalues_one_without(object))
