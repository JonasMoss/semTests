#!/usr/bin/env Rscript
# Independent fits and complete covariance construction in both libraries.
# Continuous covariance parameters are explicitly converted from magmaan's
# N-divisor fit to lavaan's (N-1)-divisor LS fit before comparison.
stopifnot(requireNamespace("magmaanlab",quietly=TRUE))
pkgload::load_all(".",quiet=TRUE)
rows<-list()
check<-function(case,quantity,actual,expected,tolerance) {
  stopifnot(length(actual)==length(expected),all(is.finite(actual)),all(is.finite(expected)))
  gap<-max(abs(actual-expected))
  rows[[length(rows)+1L]]<<-data.frame(case,quantity,gap,tolerance)
  if(gap>tolerance) stop(case," ",quantity," gap ",gap)
}
# The complete parameter sandwich that the nested score and criterion tests
# are built from, in lavaan's positional coef() order.
sandwich_vcov<-function(fit,weight_effect='estimated') {
  parts<-parameter_inference_parts(fit,'observed','centered',weight_effect,1e-4)
  inverse<-inference_solve(parts$curvature)
  B<-model_parameter_basis(fit)
  B%*%inverse%*%parts$score_covariance%*%inverse%*%t(B)
}
model<-'f =~ x1+a*x2+b*x3+x4'
dat<-lavaan::HolzingerSwineford1939[,paste0('x',1:4)]
for(est in c('GLS','WLS','ULS')) {
  fit<-lavaan::cfa(model,dat,estimator=est)
  # Each library builds its own WLS (ADF) weight; magmaan refuses a supplied
  # W for estimated-weight inference because it has no influence recipe.
  mag<-magmaanlab::fit_model(magmaanlab::model_spec(model),dat,estimator=est)
  ref<-magmaanlab::magmaan_core$infer_casewise_influence_ij_fit(mag,as.matrix(dat),NULL)
  pt<-mag$partable[mag$partable$free>0,];pt<-pt[order(pt$free),]
  scale<-ifelse(pt$op=='~~',nrow(dat)/(nrow(dat)-1),1)
  V<-crossprod(ref$influence)*tcrossprod(scale)
  Vf<-crossprod(ref$influence_naive)*tcrossprod(scale)
  check(est,'complete covariance',sandwich_vcov(fit),V,2e-6)
  check(est,'fixed covariance',sandwich_vcov(fit,'fixed'),Vf,2e-6)
}
for(categories in c(2L,3L)) {
  d<-dat
  d[]<-lapply(d,function(x) as.integer(cut(x,unique(quantile(x,seq(0,1,length.out=categories+1))),include.lowest=TRUE)))
  fit<-lavaan::cfa(model,d,ordered=names(d))
  mag<-magmaanlab::fit_model(magmaanlab::model_spec(model),d,estimator='DWLS',ordered=names(d))
  ref<-magmaanlab::magmaan_core$robust_ordinal_ij(mag,mag$ordinal_stats)
  # Both libraries use the outer-product first-stage bread, so every category
  # count is compared without substituting either side's influence.
  check(paste0('ordinal-',categories),'independent complete covariance',
        sandwich_vcov(fit),ref$vcov,2e-5)
  null<-paste(model,'a==b',sep='\n')
  fit0<-lavaan::cfa(null,d,ordered=names(d))
  mag0<-magmaanlab::fit_model(magmaanlab::model_spec(null),d,estimator='DWLS',ordered=names(d))
  # magmaan's ordinal LS score tests default to observed sensitivity
  # (bread = "observed", magmaan 189fa24a); semTests calls it curvature.
  # Compare each setting with its counterpart, covering magmaan's default and
  # the expected-curvature form.
  for(curvature in c('expected','observed')) {
    case<-paste0('ordinal-',categories,' ',curvature)
    ref0<-magmaanlab::score_tests_robust(mag0,estimated_weight=TRUE,bread=curvature)
    out<-score_test_nested(fit0,fit,curvature=curvature,weight='score_covariance')
    check(case,'nested complete studentized score',out$statistic,ref0$mi.scaled,2e-4)
    # magmaan's mi and scaling factor weight the score by expected information
    # whatever the bread; semTests' weight = 'curvature' follows the curvature.
    weighted<-score_test_nested(fit0,fit,curvature=curvature,weight='expected_curvature')
    # The magmaan score uses N; lavaan's LS criterion uses N-1.
    check(case,'nested quadratic',weighted$statistic,ref0$mi*(nrow(d)-1)/nrow(d),2e-4)
    check(case,'nested spectrum',weighted$eigenvalues,ref0$scaling.factor*(nrow(d)-1)/nrow(d),2e-4)
  }
}
out<-do.call(rbind,rows)
dir.create('workspace',showWarnings=FALSE)
write.csv(out,'workspace/magmaan-criterion.csv',row.names=FALSE)
files<-list.files(system.file(package='magmaanlab'),pattern='\\.(rdb|so)$',recursive=TRUE,full.names=TRUE)
capture.output(sessionInfo(),tools::md5sum(files),file='workspace/magmaan-criterion-session.txt')
print(out,row.names=FALSE)
