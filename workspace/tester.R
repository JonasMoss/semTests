model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50
object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")


pvalues(object, trad = NULL, eba = NULL, eba_half = c(2, 3))
