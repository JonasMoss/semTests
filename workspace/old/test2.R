#testing the non-nested restricted inaccuracy

p = 5
model <- paste("F=~", paste0("x", 1:p, collapse = "+"))
model <- paste(model, paste("\n G=~", paste0("y", 1:p, collapse = "+"), "; F~~start(0.5)*G"))

N= 10^5
p1 = .3
reps =100

one_rep <- function(seed){
  set.seed(seed)
  big.sample <- data.frame(covsim::rIG(N, sigma.hat, skewness = rep(2, 2*p),
                                       excesskurtosis=rep(20, 2*p))[[1]])
  #add groups
  big.sample$school=sample(c("a", "b"), size=N, replace=T, prob=c(p1, 1-p1))

  m0 <- cfa(model, data.frame(big.sample), group="school", estimator="MLM", group.equal="loadings")

  ugamma <- ugamma_non_nested(m0)
  us <- as.numeric(sum(diag(ugamma))/fitMeasures(m0, "df"))
  lav <- as.numeric(fitmeasures(m0,"chisq.scaling.factor"))
  c(lav, us)
}

res <- parallel::mclapply(1:reps, one_rep)
res.df <- do.call(rbind, res)
res.df <- data.frame(lav=res.df[, 1], us = res.df[, 2])

# not a good fit.
qplot(res.df$lav, res.df$us)+geom_abline()

# find problematic case

which.max(res.df$us-res.df$lav)
set.seed(157)
big.sample <- data.frame(covsim::rIG(N, sigma.hat, skewness = rep(2, 2*p),
                                     excesskurtosis=rep(20, 2*p))[[1]])
#add groups
big.sample$school=sample(c("a", "b"), size=N, replace=T, prob=c(p1, 1-p1))

m0 <- cfa(model, data.frame(big.sample), group="school", estimator="MLM", group.equal="loadings")
ugamma <- ugamma_non_nested(m0)
us <- as.numeric(sum(diag(ugamma))/fitMeasures(m0, "df"))
lav <- as.numeric(fitmeasures(m0,"chisq.scaling.factor"))
