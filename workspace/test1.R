HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

m1 <- cfa(HS.model,
          data = HolzingerSwineford1939,
          group = "school", estimator="MLM")
m0 <- cfa(HS.model,
          data = HolzingerSwineford1939,
          group = "school", estimator="MLM", group.equal="loadings")


bootstrapper(m0, m1, n_reps = 50)
semselector(m0, m1, n_reps = 50)
























#same same
lavaan:::lav_test_diff_Satorra2000(m1, m0)$scaling.factor
sum(diag(ugamma_nested(m0, m1)))/(fitmeasures(m0, "df")-fitmeasures(m1, "df"))


## TEST 2 Larger cfas

p <- 5
model <- paste("F=~", paste0("x", 1:p, collapse = "+"))
model <- paste(model, paste("\n G=~", paste0("y", 1:p, collapse = "+"), "; F~~start(0.5)*G"))
sigma.hat <- lavInspect(cfa(model, data=NULL), "sigma.hat")

# Simulate large sample
set.seed(2)
big.sample <- data.frame(covsim::rIG(10^4, sigma.hat, skewness = rep(2, 2*p),
                                     excesskurtosis=rep(20, 2*p))[[1]])
#add groups
big.sample$school=sample(c("a", "b"), size=10^4, replace=T, prob=c(.1, .9))

m0 <- cfa(model, data.frame(big.sample), group="school", estimator="MLM", group.equal="loadings")
m1 <- cfa(model, data.frame(big.sample), group="school", estimator="MLM", start=m0)

lavaan:::lav_test_diff_Satorra2000(m1, m0)$scaling.factor
sum(diag(ugamma_nested(m0, m1)))/(fitmeasures(m0, "df")-fitmeasures(m1, "df"))
#NON-NESTED TESTS:

ugamma <- ugamma_non_nested(m1)
sum(diag(ugamma))/fitMeasures(m1, "df")
fitMeasures(m1, "chisq.scaling.factor")


## not exact match.
ugamma <- ugamma_non_nested(m0)
as.numeric(sum(diag(ugamma))/fitMeasures(m0, "df"))
#/fitMeasures(m0, "df")
as.numeric(fitmeasures(m0,"chisq.scaling.factor"))


