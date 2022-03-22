rm(list=ls())

library(psych)
library(reshape2)
library(lavaan)
library(ggplot2)
library(goftest)# for the anderson-darling test. An alternative to the ks.test.


#  Here we  will use only candidate statistics already implemented in lavaan

data(bfi)# from psych package. BFI is a five-factor model (ordinal, but we treat as continuous here)
set.seed(1)
orig.sample = bfi[sample(1:nrow(bfi), size=175, replace=F), 
                  c(paste0("A", 1:5), paste0("C",1:5))] #take random sample from the large bfi-dataset
orig.sample <- orig.sample[complete.cases(orig.sample),] 

#we have a complete dataset with n=194
model ="A =~ A1+A2+A3+A4+A5; C =~ C1+C2+C3+C4+C5"


## obtain p-values from original sample. Which one to trust the most??

# 1. ML
f=sem(model, data=orig.sample, test="standard")
pml.orig <- fitmeasures(f,"pvalue") 

# 2 SB
f=sem(model, data=orig.sample, test="Satorra.Bentler")
psb.orig <- fitmeasures(f,"pvalue.scaled")

# 3 SS Scaled-and-shifted (asparouhov/muthen 2010)
f=sem(model, data=orig.sample, test="scaled.shifted")
pss.orig <- fitmeasures(f,"pvalue.scaled")

# 4 Yuan-Bentler
f=sem(model, data=orig.sample, test="Yuan.Bentler")
pyb.orig <- fitmeasures(f,"pvalue.scaled")




cat("Which of these p-values is most trusthworthy?\n ML: ", round(pml.orig, 3), "SB: ", round(psb.orig,3), "SS: ", round(pss.orig,3),
    "YB: ", round(pyb.orig,3), "\n")

#####
# transform sample a la Bollen-Stine
sigma.hat <- lavInspect(f, "sigma.hat")
sigma.sqrt <- lav_matrix_symmetric_sqrt(sigma.hat)
s.inv.sqrt <- lav_matrix_symmetric_sqrt(f@SampleStats@icov[[1]])
sample.transformed <- data.frame(as.matrix(orig.sample) %*% s.inv.sqrt %*% sigma.sqrt)
colnames(sample.transformed) <- colnames(orig.sample)
#####

#worker function. Easy to put in parallel
run.bootstrap <- function(seed, b.reps, sample.transformed){
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = b.reps, style = 3)
  pml=NULL; psb=NULL;  pss=NULL; pyb=NULL
  for(b in 1:b.reps){
    setTxtProgressBar(pb, b)
    boot.sample = sample.transformed[sample(1:nrow(sample.transformed), replace=T), ]
    f.sb = tryCatch(sem(model,boot.sample, test="Satorra.Bentler", start=f), error=function(w) { TRUE})
    f.ss = tryCatch(sem(model, boot.sample, test="scaled.shifted", start=f.sb),error=function(w) { TRUE})
    f.yb = tryCatch(sem(model,boot.sample, test="Yuan.Bentler", start=f.sb), error=function(w) { TRUE})
    
    while(isTRUE(f.sb) | isTRUE(f.ss) |  isTRUE(f.yb)  ){
      boot.sample = sample.transformed[sample(1:nrow(sample.transformed), replace=T), ]
      f.sb = tryCatch(sem(model,boot.sample, test="Satorra.Bentler", start=f), error=function(w) { TRUE})
      f.ss = tryCatch(sem(model, boot.sample, test="scaled.shifted", start=f.sb),error=function(w) { TRUE})
      f.yb = tryCatch(sem(model,boot.sample, test="Yuan.Bentler", start=f.sb), error=function(w) { TRUE})
    }
    
    pml <- c(pml,tryCatch(fitmeasures(f.sb,"pvalue"), error=function(w) { TRUE}))
    psb <- c(psb,tryCatch(fitmeasures(f.sb,"pvalue.scaled"), error=function(w) { TRUE}))
    pss <- c(pss,tryCatch(fitmeasures(f.ss,"pvalue.scaled"), error=function(w) { TRUE}))
    pyb <- c(pyb,tryCatch(fitmeasures(f.yb,"pvalue.scaled"), error=function(w) { TRUE}))
  }
  return(data.frame(pml,psb,pss,  pyb))
}

## bootstrap reps
b.reps=1000

## run bootstrap in serial 
res.df <- run.bootstrap(X=1, b.reps=b.reps, sample.transformed=sample.transformed)

colnames(res.df)<- c("ML", "SB", "SS",  "YB")
rownames(res.df) = NULL 
melted = melt(res.df)
melted$statistic <- melted$variable

#plot of the p-values
ggplot(melted, aes(value, fill=statistic))+geom_histogram()+facet_wrap(~statistic)

#Anderson-Darling distance to uniform distribution
library(goftest)
distUniform= as.vector(sapply(my.df,function(x) ad.test(x, null="punif")$statistic ))
distUniform

cat("Anderson-Darling test selects: ", colnames(my.df)[which.min(distUniform)])

# Look at the rejection rates for significance levels between 0.01 and 0.1:
x=seq(0.01, 0.1, length.out=10)
y <- t(sapply(x, function(x) colMeans(my.df<x)))
m <- melt(data.frame(cbind(x, y)), id.vars="x"); m$test <- m$variable
ggplot(m, aes(x, value, color=test))+geom_line()+geom_abline(slope=1, intercept=0)+ylab("Rejection rate")+
  xlab("Level of significance")+geom_point()
