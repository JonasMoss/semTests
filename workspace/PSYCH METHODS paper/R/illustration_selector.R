library(CompQuadForm)

source("eigenAux.R")

model <- "A =~ A1+A2+A3+A4+A5;
         C =~ C1+C2+C3+C4+C5"

rownumber <- 200
t <- psych::bfi[1:rownumber, 1:10]
f <- sem(model, t)



# transform sample for bootstrapping draws
Sigma.hat <- lavaan:::computeSigmaHat(lavmodel = f@Model)
sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat[[1]])
S.inv.sqrt <- lav_matrix_symmetric_sqrt(f@SampleStats@icov[[1]])
sample.transformed <- data.frame(as.matrix(t) %*% S.inv.sqrt %*% sigma.sqrt)
colnames(sample.transformed) <- colnames(t)

tml <- fitmeasures(f, "chisq")
pml <- fitmeasures(f, "pvalue")

DF <- fitMeasures(f, "df")
UG <- inspect(f, "UG")
lambdas <- Re(eigen(UG)$values[1:DF])

eigenps <- get.pvalues.new(tml, lambdas) # full, half and sb

p.cf <- scaled.F(tml, lambdas)
# asparouhov scaled and shifted
pss <- fitmeasures(sem(model, t, test = "scaled.shifted"), "pvalue.scaled")


## we use: pml, psb, pss, pfull, phalf, pcf

R <- 5000
set.seed(1)
pml <- NULL
psb <- NULL
pss <- NULL
pfull <- NULL
phalf <- NULL
pcf <- NULL
for (i in 1:R) {
  if (i %% 51 == 0) {
    print(i)
  }
  boot.sample <- sample.transformed[sample(1:rownumber, replace = T), ]

  f.boot <- tryCatch(sem(model, boot.sample), error = function(w) {
    TRUE
  })
  UG <- tryCatch(inspect(f.boot, "UG"), error = function(w) {
    TRUE
  })
  while (isTRUE(f.boot) || !(inspect(f.boot, "converged")) || isTRUE(UG)) {
    cat("error ")
    boot.sample <- sample.transformed[sample(1:rownumber, replace = T), ]
    f.boot <- tryCatch(sem(model, boot.sample), error = function(w) {
      TRUE
    })
    UG <- tryCatch(inspect(f.boot, "UG"), error = function(w) {
      TRUE
    })
  }


  tml.boot <- fitmeasures(f.boot, "chisq")
  pml <- c(pml, fitmeasures(f.boot, "pvalue"))

  DF <- fitMeasures(f.boot, "df")
  lambdas <- Re(eigen(UG)$values[1:DF])

  eigenps.boot <- get.pvalues.new(tml.boot, lambdas) # full, half and sb
  psb <- c(psb, unlist(eigenps.boot[3]))
  pfull <- c(pfull, unlist(eigenps.boot[1]))
  phalf <- c(phalf, unlist(eigenps.boot[2]))

  pcf <- c(pcf, scaled.F(tml.boot, lambdas))

  # asparouhov scaled and shifted
  pss <- c(pss, fitmeasures(sem(model, boot.sample, test = "scaled.shifted"), "pvalue.scaled"))
}


my.df <- data.frame(pml, psb, pss, pfull, phalf, pcf)
saveRDS(my.df, "illustration1_R5000.rds")
theme_set(theme_minimal())


colnames(my.df) <- c("ML", "SB", "SS", "EBAF", "EBA2", "CF")

ggplot(melt(my.df), aes(value)) +
  geom_histogram() +
  xlab("p-value") +
  facet_wrap(~variable, nrow = 3)

# KS distance
distUniform <- as.vector(sapply(my.df, function(x) ks.test(x, y = "punif")$statistic))
cat("KS test selects: ", colnames(my.df)[which.min(distUniform)])


# AD
distUniform <- as.vector(sapply(my.df, function(x) ad.test(x, null = "punif")$statistic))
distUniform
cat("Anderson-Darling test selects: ", colnames(my.df)[which.min(distUniform)])
