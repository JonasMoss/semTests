model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5;
          E =~ E1+E2+E3+E4+E5;
          N =~ N1+N2+N3+N4+N5;
          O =~ O1+O2+O3+O4+O5;"
obj <- lavaan::sem(model, na.omit(psych::bfi), estimator = "MLM")

obj <- lavaan::sem(model, psych::bfi, estimator = "MLM")
x <- psych::bfi


now <- Sys.time()
obj <- lavaan::cfa(model, data = x, estimator = "ML", std.lv = TRUE)
#gamma <- lavaan::lavInspect(obj, what ="gamma")
lavaan::lavInspect(obj, what ="UGamma")
later <- Sys.time()
later - now



now <- Sys.time()
obj <- lavaan::cfa(model, data = x, estimator = "ML", std.lv = TRUE)
object@Options$gamma.unbiased <- TRUE
u <- lavaan::lavInspect(obj, what ="UfromUGamma")
gamma <- lavaan::lavInspect(obj, what ="gamma")
res <- u %*% gamma
later <- Sys.time()
later - now


object@Options$gamma.unbiased <- FALSE
now <- Sys.time()
obj <- lavaan::cfa(model, data = x, estimator = "ML", std.lv = TRUE)
res <- lavaan::lavInspect(obj, what ="UGamma")
later <- Sys.time()
later - now


now <- Sys.time()
obj <- lavaan::cfa(model, data = x, estimator = "MLM", std.lv = TRUE)
#gamma <- lavaan::lavInspect(obj, what ="gamma")
later <- Sys.time()
later - now

now <- Sys.time()
uu <- lavaan::lavInspect(obj, what ="UfromUGamma")
later <- Sys.time()
later - now

now <- Sys.time()
u <- lavaan::lavInspect(obj, what ="gamma")
later <- Sys.time()
later - now

object@Options$gamma.unbiased <- unbiased


now <- Sys.time()
ug <- lavaan::lavInspect(obj, what ="UGamma")
later <- Sys.time()
later - now



str(obj)

hyp <- obj@SampleStats@NACOV
o


sum(abs(hyp[[1]] / nrow(na.omit(x)) - gamma[[1]]))

now <- Sys.time()
lavaan::lavInspect(obj, what ="U")
later <- Sys.time()

now <- Sys.time()
gamma <- lavaan::lavInspect(obj, what ="gamma")
later <- Sys.time()
later - now
