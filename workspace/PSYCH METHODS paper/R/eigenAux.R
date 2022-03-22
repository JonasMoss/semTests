## help functions for eigenvalue Testing of SEM
library(lavaan)
require(CompQuadForm) ## The "Imhof"-method is described as exact in the author's comparison-paper.
require(MASS)

### MODELS
bollen.model <- '
# measurement model
ind60 =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
y4 ~~ y8
y6 ~~ y8
'
bollen.population.model <- '
# latent variable definitions
ind60 =~ x1 + 0.8*x2 + 0.3*x3
dem60 =~ y1 + 0.8*y2 + y3 + y4
dem65 =~ y5 + y6 + 0.8*y7 + y8
# regressions
dem60 ~ 0.3*ind60
dem65 ~ 0.2*ind60 + 0.8*dem60
# residual (co)variances
y1 ~~ .3*y5
y2 ~~ .4* y4 + .1*y6
y3 ~~ .3*y7
y4 ~~ .3*y8
y6 ~~ .3*y8
'
bollen.restricted <- '
# measurement model
ind60 =~ x1 + a*x2 + x3
dem60 =~ y1 + a*y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ c*y5
y2 ~~ y4 + y6
y3 ~~ c*y7
y4 ~~ c*y8
y6 ~~ c*y8
'

### IN PAPER: 13 df.
restricted.testing <- '
# measurement model
ind60 =~ x1 + a*x2 + x3
dem60 =~ y1 + a*y2 + y3 + d*y4
dem65 =~ y5 + y6 + y7 + d*y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ c*y5
y2 ~~ y4 + y6
y3 ~~ c*y7
y4 ~~ c*y8
y6 ~~ c*y8
y1~~ e*y1
y2~~ e*y2
y3~~ e*y3
y4~~ e*y4
y5~~ e*y5
y6~~ e*y6
y7~~ e*y7
y8~~ e*y8

dem60 ~~ b *dem60
dem65 ~~ b *dem65

'




bollen.model.lackstworesvars <- '
# measurement model
ind60 =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ y5
y2 ~~  y6
y3 ~~ y7
y4 ~~ y8

'

bollen.model.lackmoderate <- '
# measurement model
ind60 =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
'

bollen.model.lacksevere <- '
# measurement model
ind60 =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
'
sigma0 <- inspect(sem(bollen.population.model,  data=NULL), "sigma.hat")

### DWLS functions.


