library("semselector")
library("future.apply")
plan(multisession)

set.seed(313)
model <- "A =~ A1+A2+A3+A4+A5"
object <- lavaan::sem(model, psych::bfi[1:150, ], estimator = "ML")
pvalues(object)



model <- "A =~ A1+A2+A3+A4+A5"
model <- "A =~ A1+A2+A3+A4+A5;
          O =~ O1+O2+O3+O4+O5;
          E =~ E1+E2+E3+E4+E5;
          N =~ N1+N2+N3+N4+N5;
          C =~ C1+C2+C3+C4+C5"

lavaan::sem(model, psych::bfi, estimator = "ML")


model <- "A =~ A1+A2+A3+A4+A5+O1+O2+O3+O4+O5+E1+E2+E3+E4+E5+N1+N2+N3+N4+N5+C1+C2+C3+C4+C5"


n <- 200
object1 <- lavaan::sem(model, psych::bfi, estimator = "ML")
object2 <- lavaan::sem(model, psych::bfi, estimator = "MLM")



microbenchmark::microbenchmark(
  lavaan::sem(model, psych::bfi, estimator = "ML"),
  lavaan::sem(model, psych::bfi, estimator = "MLM")
  , times = 10)

microbenchmark::microbenchmark(
  pvalues_one(object1),
  pvalues_one(object2)
  , times = 10
)

microbenchmark::microbenchmark(
  pvalues_one(lavaan::sem(model, psych::bfi, estimator = "ML")),
  pvalues_one(lavaan::sem(model, psych::bfi, estimator = "MLM"))
  , times = 10)
