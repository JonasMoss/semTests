test_that("observed saturated curvature includes fitted-mean residual terms", {
  S <- matrix(c(2,.4,.4,1.5),2)
  mu <- c(.7,-.4)
  X <- rbind(c(-1,.2),c(2,1),c(.3,2.1))
  geo <- gaussian_moment_geometry(S,TRUE)
  par <- c(mu,S[geo$index])
  grad <- function(x) {
    V <- matrix(0,2,2); V[geo$index] <- x[-(1:2)]
    V[upper.tri(V)] <- t(V)[upper.tri(V)]
    colSums(gaussian_score_rows(sweep(X,2,x[1:2],"-"),gaussian_moment_geometry(V,TRUE)))
  }
  H <- criterion_hessian(grad,par,diag(5),1e-5)$matrix
  expect_equal(gaussian_observed_information(sweep(X,2,mu,"-"),geo),H,tolerance=1e-8)
})

test_that("profile geometry is invariant to nuisance and alternative coordinates", {
  set.seed(51)
  Q <- matrix(rnorm(25),5)
  H <- crossprod(Q)+diag(5)
  J <- crossprod(matrix(rnorm(35),7,5))+diag(5)
  K <- diag(5)[,1:3]
  s <- c(0,0,0,.3,-.8)
  a <- profile_score_test_matrices(s,H,J,K,keep=TRUE)
  T <- diag(c(.5,2,3,1,4))
  b <- profile_score_test_matrices(drop(crossprod(T,s)),crossprod(T,H%*%T),
    crossprod(T,J%*%T),solve(T,K%*%diag(c(2,3,.5))))
  expect_equal(a$statistic,b$statistic,tolerance=1e-10)
  expect_equal(a$eigenvalues,b$eigenvalues,tolerance=1e-10)
  expect_lt(max(abs(crossprod(K,H%*%a$ingredients$projection))),1e-12)
  expect_error(profile_score_test_matrices(s,H,J,cbind(K,K[,1])),"rank deficient")
})

test_that("FIML pattern projection reduces to complete ML and handles missingness", {
  d <- lavaan::HolzingerSwineford1939
  model <- 'f =~ x1+x2+x3+x4'
  f <- lavaan::cfa(model,d,meanstructure=TRUE)
  fm <- lavaan::cfa(model,d,missing="ml")
  a <- score_test(f); b <- score_test(fm)
  expect_equal(a$statistic,b$statistic,tolerance=1e-5)
  expect_equal(a$eigenvalues,b$eigenvalues,tolerance=1e-5)
  set.seed(112)
  d$x2[sample(nrow(d),60)] <- NA
  d$x4[sample(nrow(d),50)] <- NA
  fm <- lavaan::cfa(model,d,missing="ml")
  patterns <- gaussian_pattern_parts(fm)
  all <- pattern_score_products(patterns,diag(nrow(model_delta_matrix(fm))),"centered")
  G <- score_test(fm,keep=TRUE,score_covariance ="centered")$ingredients$projection
  projected <- pattern_score_products(patterns,G,"centered",tile_rows=1L)
  expect_equal(projected$score,drop(crossprod(G,all$score)),tolerance=1e-9)
  expect_equal(projected$score_covariance,crossprod(G,all$score_covariance%*%G),tolerance=1e-9)
  expect_lt(score_test(fm)$diagnostics$nuisance_stationarity_norm,1e-4)
  expect_true(is.finite(score_test(fm,curvature ="observed")$p_value))
})

test_that("complete ML sandwich uses the actual model Hessian and scores", {
  d <- lavaan::HolzingerSwineford1939
  d$x1 <- d$x1-mean(d$x1)+.3
  f <- lavaan::cfa('f =~ x1+x2+x3+x4\nx1 ~ 0*1',d,meanstructure=TRUE)
  parts <- parameter_inference_parts(f)
  B <- model_parameter_basis(f)
  # lavaan's numerical objective Hessian is an independent criterion route.
  ref <- get("lav_model_hessian",asNamespace("lavaan"))(f@Model,f@SampleStats,f@Data,f@Options,f@Cache)
  expect_equal(unname(parts$curvature/f@SampleStats@ntotal),unname(crossprod(B,ref%*%B)),tolerance=2e-5)
  half <- parameter_inference_parts(f,derivative_step=5e-5)
  expect_equal(parts$curvature,half$curvature,tolerance=2e-6)
  expect_gt(max(abs(parts$curvature-parts$expected)),.01)
  R <- matrix(0,1,length(lavaan::coef(f))); R[1,1] <- 1
  w <- wald_test(f,R,vcov="sandwich",keep=TRUE)
  expected <- inference_solve(parts$curvature)%*%parts$score_covariance%*%inference_solve(parts$curvature)
  expect_equal(as.numeric(w$ingredients$vcov),as.numeric(B%*%expected%*%t(B)),tolerance=1e-9)
})

test_that("nested score accepts an unoptimized affine template and checks embedding", {
  d <- lavaan::HolzingerSwineford1939
  model <- 'f =~ x1+a*x2+b*x3+x4'
  f0 <- lavaan::cfa(paste(model,'a==b',sep='\n'),d)
  f1 <- lavaan::cfa(model,d)
  template <- lavaan::cfa(model,d,do.fit=FALSE)
  a <- score_test_nested(f0,f1,keep=TRUE)
  b <- score_test_nested(f0,template)
  expect_equal(a$statistic,b$statistic,tolerance=1e-10)
  expect_equal(a$eigenvalues,b$eigenvalues,tolerance=1e-10)
  expect_equal(a$df,1L)
  expect_lt(a$diagnostics$nuisance_stationarity_norm,1e-5)
  expect_identical(score_test_nested(f0,f1,metric="sandwich")$reference,"chisq")
  lr <- criterion_test_nested(f0,f1)
  expect_equal(lr$statistic,as.numeric(lavaan::fitMeasures(f0,"chisq")-lavaan::fitMeasures(f1,"chisq")),tolerance=1e-9)
  expect_equal(lr$eigenvalues,a$eigenvalues,tolerance=1e-10)
  expect_error(score_test_nested(f1,f0),"nested|constraints|directions")
  expect_error(score_test_nested(f0,f0),"directions")
})

test_that("continuous complete weight influences agree with empirical contamination", {
  d <- lavaan::HolzingerSwineford1939
  for (est in c("GLS","WLS","DWLS","ULS")) {
    f <- lavaan::cfa('f =~ x1+x2+x3+x4',d,estimator=est,ordered=FALSE,meanstructure=TRUE)
    first <- continuous_first_stage(f,1,"estimated")
    X <- f@Data@X[[1]]
    n <- nrow(X); p <- ncol(X)
    ij <- which(lower.tri(matrix(0,p,p),diag=TRUE),arr.ind=TRUE)
    fitted <- inference_implied(f)[[1]]
    # Recompute the entire empirical functional, independently of analytic IJ.
    functional <- function(prob) {
      mu <- colSums(X*prob)
      xc <- sweep(X,2,mu,"-")
      S <- crossprod(xc,xc*prob)
      z <- cbind(xc,sweep(xc[,ij[,1]]*xc[,ij[,2]],2,S[ij],"-"))
      Gamma <- crossprod(z,z*prob)
      nt <- gaussian_moment_geometry(S*n/(n-1),TRUE)$information
      V <- switch(est,GLS=nt,WLS=solve(Gamma),DWLS=diag(1/diag(solve(nt))),ULS=diag(ncol(z)))
      moment <- c(mu,S[ij]*n/(n-1))
      drop(V%*%(moment-fitted))
    }
    for (i in c(1L,23L)) {
      direction <- rep(-1/n,n);direction[i] <- direction[i]+1
      h <- 1e-6
      fd <- (functional(rep(1/n,n)+h*direction)-functional(rep(1/n,n)-h*direction))/(2*h)
      actual <- (first$influence%*%first$weight+first$correction(first$observed-fitted))[i,]
      expect_equal(as.numeric(actual),fd,tolerance=2e-6)
    }
    expect_equal(unname(first$correction(rep(0,length(fitted)))),matrix(0,n,length(fitted)))
    parts <- parameter_inference_parts(f)
    expect_lt(max(abs(parts$score))/n,1e-4)
    expect_true(all(is.finite(parts$score_covariance)))
  }
})

test_that("ordinal numerical IJ converges and estimated weights matter off null", {
  d <- lavaan::HolzingerSwineford1939[,paste0("x",1:4)]
  d[] <- lapply(d,function(x) as.integer(x>median(x)))
  f <- lavaan::cfa('f =~ x1+x2+x3+x4',d,ordered=names(d))
  first <- ordinal_first_stage(f,1,"estimated",1e-4)
  half <- ordinal_first_stage(f,1,"estimated",5e-5)
  r <- first$observed-inference_implied(f)[[1]]
  expect_equal(first$influence,half$influence,tolerance=1e-5)
  expect_equal(first$correction(r),half$correction(r),tolerance=1e-5)
  expect_gt(max(abs(first$correction(r))),1e-4)
  expect_equal(first$correction(0*r),matrix(0,nrow(d),length(r)))
  R <- matrix(0,1,length(lavaan::coef(f)));R[1,1]<-1
  out <- wald_test(f,R,vcov="sandwich",keep=TRUE)
  fixed <- wald_test(f,R,vcov="sandwich",weight_effect="fixed")
  expect_true(is.finite(out$p_value))
  expect_gt(abs(out$statistic-fixed$statistic),1e-5)
  expect_error(wald_test(f,R,vcov="sandwich",score_covariance ="uncentered"),"centered")
})

test_that("multi-category first-stage derivatives are stable and retain native weights", {
  d <- lavaan::HolzingerSwineford1939[,paste0("x",1:4)]
  d[] <- lapply(d,function(x) as.integer(cut(x,quantile(x,c(0,1/3,2/3,1)),include.lowest=TRUE)))
  f <- lavaan::cfa('f =~ x1+x2+x3+x4',d,ordered=names(d))
  a <- ordinal_first_stage(f,1,"estimated",1e-4)
  b <- ordinal_first_stage(f,1,"estimated",5e-5)
  r <- a$observed-inference_implied(f)[[1]]
  expect_equal(a$influence,b$influence,tolerance=1e-7)
  expect_equal(a$correction(r),b$correction(r),tolerance=1e-7)
  # Marginal quantile influences supply an independent analytic check.
  thresholds <- a$observed[1:8]
  q <- 0L
  for(j in 1:4) for(k in 1:2) {
    q <- q+1L
    expected <- (as.numeric(f@Data@X[[1]][,j]<=k)-stats::pnorm(thresholds[q]))/stats::dnorm(thresholds[q])
    expect_equal(a$influence[,q],expected,tolerance=1e-7)
  }
  # The influence uses the outer-product bread, which is lavaan's categorical
  # Gamma convention, so it reproduces NACOV rather than merely approximating
  # it. Under latent bivariate normality the pair model is correctly specified,
  # the information identity holds, and this agrees asymptotically with the
  # exact-Jacobian bread. The two differ in finite samples only once a table
  # exceeds 2x2, where (a-1)+(b-1)+1 parameters stop saturating its ab-1 cells.
  expect_equal(crossprod(a$influence)/nrow(d), f@SampleStats@NACOV[[1]],
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("unequal-group LS curvature follows the fitted N-minus-one criterion", {
  f <- lavaan::cfa('f =~ x1+x2+x3+x4',lavaan::HolzingerSwineford1939,
                  estimator="GLS",group="school",group.equal="loadings",meanstructure=TRUE)
  a <- parameter_inference_parts(f)
  H <- get("lav_model_hessian",asNamespace("lavaan"))(f@Model,f@SampleStats,f@Data,f@Options,f@Cache)
  B <- model_parameter_basis(f)
  expect_equal(unname(a$curvature/a$nobs),unname(crossprod(B,H%*%B)),tolerance=2e-5)
  expect_lt(max(abs(a$score))/a$nobs,1e-5)
})

test_that("FIML sufficient-statistic gradients match casewise scores and curvature", {
  set.seed(880)
  d <- lavaan::HolzingerSwineford1939
  d$x1 <- d$x1-mean(d$x1)+.2
  d$x2[sample(nrow(d),60)] <- NA
  d$x4[sample(nrow(d),50)] <- NA
  f <- lavaan::cfa('f =~ x1+x2+x3+x4\nx1 ~ 0*1',d,missing="ml")
  B <- model_parameter_basis(f)
  theta <- as.numeric(lavaan::coef(f))
  records <- gaussian_pattern_parts(f)$records
  perturbed <- inference_model_at(f,theta+c(.001,rep(0,length(theta)-1)))
  expect_equal(gaussian_parameter_gradient(perturbed,B,records),
               gaussian_parameter_evaluation(perturbed,B)$score,tolerance=1e-9)
  parts <- parameter_inference_parts(f)
  H <- get("lav_model_hessian",asNamespace("lavaan"))(f@Model,f@SampleStats,f@Data,f@Options,f@Cache)
  expect_equal(unname(parts$curvature/parts$nobs),unname(crossprod(B,H%*%B)),tolerance=2e-5)
})

test_that("unoptimized templates use affine offsets rather than starting estimates", {
  d <- lavaan::HolzingerSwineford1939
  model <- 'f =~ x1+a*x2+b*x3+c*x4\na==2*b+.1'
  f0 <- lavaan::cfa(paste(model,'c==b',sep='\n'),d,estimator="MLM")
  template <- lavaan::cfa(model,d,do.fit=FALSE)
  f1 <- lavaan::cfa(model,d)
  expect_gt(max(abs(template@Model@ceq.function(lavaan::coef(template)))),.01)
  a <- score_test_nested(f0,template)
  b <- score_test_nested(f0,f1)
  expect_equal(a$statistic,b$statistic,tolerance=1e-10)
  expect_equal(a$eigenvalues,b$eigenvalues,tolerance=1e-10)
})
