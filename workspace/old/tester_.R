library("semselector")
set.seed(313)
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 200
object <- lavaan::sem(model, psych::bfi[1:n, 1:10])
pvalues(object)

selector <- semselector(object, n_reps = 500)
attr(selector, "bollen-stine")





m0 = lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "ML")

m1 = lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")

sapply(seq(m0@Options), function(i) m0@Options[[i]] == m1@Options[[i]])
