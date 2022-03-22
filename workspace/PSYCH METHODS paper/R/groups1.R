HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

f0 <- cfa(HS.model, 
           data = HolzingerSwineford1939, 
            estimator="MLM")
str(lavInspect(f0,"Gamma"))

lavInspect(f0,"Gamma")
str(lavInspect(f0,"UGamma"))

fit <- cfa(HS.model, 
           data = HolzingerSwineford1939, 
           group = "school", estimator="MLM")
str(lavInspect(fit, "Gamma"))
str(lavInspect(fit, "UGamma"))

fit <- cfa(HS.model, 
           data = HolzingerSwineford1939, 
           group = "school", estimator="MLM")

lavTestLRT(f0, fit)
anova(f0, fit)

anova(cfa(HS.model, 
    data = HolzingerSwineford1939), cfa(HS.model, 
                                        data = HolzingerSwineford1939, group="school"))

## intercept
fit1 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school", estimator="MLM")

# weak invariance
fit2 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school", estimator="MLM",
            group.equal = "loadings")

lavTestLRT(fit1, fit2)
