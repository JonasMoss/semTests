# the slowness of pvalues()
library("lavaan")
rm(list=ls())
model <- paste0("F=~", paste0("x",1:40, collapse="+"))
set.seed(1)
data <- lavaan::simulateData(model)
fit <- lavaan::cfa(model, data, estimator="MLM", sample.nobs=1000)

system.time(cfa(model, data))#.3 sec
system.time(cfa(model, data, estimator="MLM",std.lv=T))#1.5 seec

system.time(cfa(model, data, test="scaled.shifted",std.lv=T))#2.3 sec
system.time(semTests::pvalues(fit, tests="SB"))#3.5 sec NO Warning when ML?
system.time(semTests::pvalues(fit, tests=c("SB","SB_RLS")))#8.3!!
system.time(ebas <- semTests::pvalues(fit))#24 why so slow?


system.time(semTests::pvalues(fit, tests=c("SB")))#8.3!!
system.time(semTests::pvalues(fit, tests=c("SB_UG", "SB")))#8.3!!

system.time(pvalues(fit, tests=c("SB")))#8.3!!
system.time(pvalues(fit, tests=c("SB_UG", "SB")))#8.3!!


system.time(semTests::pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = NULL, pols = NULL, unbiased = 1, chis = c("ml")))#8.3!!
system.time(semTests::pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = NULL, pols = NULL, unbiased = 2, chis = c("ml")))#8.3!!
system.time(semTests::pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = NULL, pols = NULL, unbiased = 3, chis = c("ml")))#8.3!!

system.time(pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = NULL, pols = NULL, unbiased = 1, chis = c("ml")))#8.3!!
system.time(pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = NULL, pols = NULL, unbiased = 2, chis = c("ml")))#8.3!!
system.time(pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = NULL, pols = NULL, unbiased = 3, chis = c("ml")))#8.3!!


system.time(semTests::pvalues(fit, tests=c("SB_RLS")))#8.3!!




system.time(semTests::pvalues(fit, tests=c("SB_RLS_UG")))#8.3!!
system.time(semTests::pvalues(fit, tests=c("SB_UG","SB_RLS_UG")))#8.3!!

system.time(semTests::pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = c(2, 4, 6), pols = 2, unbiased = 3, chis = c("ml", "rls")))#8.3!!
system.time(pvalues(fit, tests = NULL, trad = c("sb", "ss"), peba = c(2, 4, 6), pols = 2, unbiased = 3, chis = c("ml", "rls")))#8.3!!



#imhof is not the culprit
#speed and imhof 740 eigenvalues
UG <- lavInspect(fit, "UG")
lambdas <- Re(eigen(UG)$values)[1:740]# 1 sec

#compquadform is very fast
system.time(res1 <- CompQuadForm::imhof(sum(lambdas), lambdas))

#from lavaan fork
lav_fmg_imhof <- \(x, lambda) {
  theta <- \(u, x, lambda) 0.5 * (sum(atan(lambda * u)) - x * u)
  rho <- \(u, lambda) exp(1 / 4 * sum(log(1 + lambda^2 * u^2)))
  integrand <- Vectorize(\(u) {
    sin(theta(u, x, lambda)) / (u * rho(u, lambda))
  })
  z <- integrate(integrand, lower = 0, upper = Inf)$value
  0.5 + z / pi
}
system.time(res2 <- lav_fmg_imhof(sum(lambdas), lambdas))
