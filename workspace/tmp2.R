x1 <- lav_ugamma_nested_2000(m0, m1)
x2 <- lav_ugamma_nested_2001(m0, m1)

x = sort(Re(eigen(x1)$values))
y = sort(Re(eigen(x2)$values))

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

pvalues_nested(m0, m1, tests = "SB_UG")
pvalues_nested(m0, m1, tests = "SB_UG_RLS")
pvalues_nested(m0, m1, tests = "SB_RLS")
pvalues_nested(m0, m1, tests = "SB")

sum(abs(ugamma_nested(m0, m1, "2000", 2)[[1]] - ugamma_nested(m0, m1, "2000", 1)[[1]]))

-sort(-Re(eigen(x2)$values))[1:6]
m1
min(Re(eigen(x2)$values))


make_chisqs(c("rls"), m0)

lavaan::lavTest(object, test = "browne.residual.nt")


m0@Options$gamma.unbiased <- FALSE
unbiased <- lavaan:::lav_object_inspect_UGamma(m0)
m0@Options$gamma.unbiased <- TRUE
biased <- lavaan:::lav_object_inspect_UGamma(m0)

m0@Options$gamma.unbiased <- TRUE
out <- lavaan:::lav_test_satorra_bentler(lavobject     = m0,
                                method        = "original",
                                return.ugamma = TRUE)
biased <- out$UGamma

out <- lavaan:::lav_test_satorra_bentler(lavobject     = m0,
                                method        = "original",
                                return.u      = TRUE)
u = out$UfromUGamma



lavoptions = lavaan::lavInspect(object, "options")
lavdata = object@Data
lavoptions$gamma.unbiased <- TRUE

#| Extract gamma from lavaan object.
get_gamma <- \(object, unbiased = FALSE) {
  lavoptions = lavaan::lavInspect(object, "options")
  lavdata = object@Data
  lavoptions$gamma.unbiased <- unbiased
  gamma_list <- lavaan:::lav_samplestats_from_data(lavdata, lavoptions)@NACOV
  lavaan:::lav_matrix_bdiag(gamma_list)
}

sum(abs(get_gamma(object, TRUE) - get_gamma(object, FALSE)))

u <- lavaan:::lav_matrix_bdiag(lavaan::lavInspect(object, "gamma"))
sum(u - ugamma)
