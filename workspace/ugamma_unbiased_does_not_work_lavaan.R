model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

pvalues_one(object)

pvalues(object)

list(ptrad = c("pstd", "psb", "pss"), peba = c(2:4))





model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")


u <- lavaan::lavInspect(object, "UGamma")
object@Options$gamma.unbiased = TRUE
gamma_unbiased <- lavaan:::lav_object_gamma(object)[[1]]
ugamma_unbiased <- lavaan:::lav_object_inspect_UGamma(object)
ugamma2_unbiased <- u %*% gamma_unbiased

object@Options$gamma.unbiased = FALSE
gamma_biased <- lavaan:::lav_object_gamma(object)[[1]]
ugamma_biased <- lavaan:::lav_object_inspect_UGamma(object)
ugamma2_biased <- u %*% gamma_biased

c(gamma_unbiased - gamma_biased)
c(ugamma_unbiased[[1]] - ugamma_biased[[1]])
c(ugamma2_unbiased - ugamma2_biased)

lavaan::lavInspect(object, "Ugamma")









sel <- semselector(object, n_reps = 2, distances = "anderson-darling", unbiased = 3)
attr(sel, "pvalues")
pvalues(object, unbiased = 3)
