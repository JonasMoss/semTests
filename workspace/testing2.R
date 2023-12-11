library("semTests")

## One model, no restrictions
set.seed(313)
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

ugamma_1 <- ugamma_non_nested(object, unbiased = TRUE)
ugamma_2 <- ugamma_non_nested(object, unbiased = FALSE)
sum(abs(ugamma_1 - ugamma_2))


object@Options$gamma.unbiased = FALSE
gamma <- lavaan:::lav_object_gamma(object)
ugamma <- lavaan:::lav_object_inspect_UGamma(object)
u <- lavaan:::lav_object_inspect_UfromUGamma(object)
u2 <- lavaan::lavInspect(object, "U")

sum(abs(u %*% gamma[[1]] - ugamma))
