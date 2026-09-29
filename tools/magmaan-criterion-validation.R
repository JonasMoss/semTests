#!/usr/bin/env Rscript
# Independent fits and complete covariance construction in both libraries.
# Continuous covariance parameters are explicitly converted from magmaan's
# N-divisor fit to lavaan's (N-1)-divisor LS fit before comparison.
stopifnot(requireNamespace("magmaan",quietly=TRUE))
pkgload::load_all(".",quiet=TRUE)
rows<-list()
diagnostics<-list()
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
  W<-if(est=='WLS') as.matrix(lavaan::lavInspect(fit,'wls.v')) else NULL
  mag<-magmaan::magmaan(magmaan::model_spec(model),dat,estimator=est,W=W)
  ref<-magmaan::magmaan_core$infer_casewise_influence_ij_fit(mag,as.matrix(dat),W)
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
  mag<-magmaan::magmaan(magmaan::model_spec(model),d,estimator='DWLS',ordered=names(d))
  ref<-magmaan::magmaan_core$robust_ordinal_ij(mag,mag$ordinal_stats)
  vcov<-sandwich_vcov(fit)
  if(categories>2L) {
    # magmaan stores an OPG first-stage influence. semTests differentiates the
    # actual two-step estimator. Record the native gap, then align ONLY this
    # first-stage influence to isolate the independently coded weight calculus.
    diagnostics[[length(diagnostics)+1L]] <- data.frame(categories,
      native_covariance_gap=max(abs(vcov-ref$vcov)))
    first<-ordinal_first_stage(fit,1,"estimated",1e-4)
    mag$ordinal_stats$moment_influence[[1]]<-first$influence
    ref<-magmaan::magmaan_core$robust_ordinal_ij(mag,mag$ordinal_stats)
  }
  label<-if(categories==2L) "independent complete covariance" else "complete covariance with aligned first-stage influence"
  check(paste0('ordinal-',categories),label,vcov,ref$vcov,2e-5)
  null<-paste(model,'a==b',sep='\n')
  fit0<-lavaan::cfa(null,d,ordered=names(d))
  mag0<-magmaan::magmaan(magmaan::model_spec(null),d,estimator='DWLS',ordered=names(d))
  if(categories>2L) {
    first0<-ordinal_first_stage(fit0,1,"estimated",1e-4)
    mag0$ordinal_stats$moment_influence[[1]]<-first0$influence
  }
  ref0<-magmaan::score_tests_robust(mag0,estimated_weight=TRUE)
  out<-score_test_nested(fit0,fit,curvature ='expected',metric='sandwich')
  check(paste0('ordinal-',categories),'nested complete studentized score',out$statistic,ref0$mi.scaled,2e-4)
  weighted<-score_test_nested(fit0,fit,curvature ='expected')
  # The magmaan score uses N; lavaan's LS criterion uses N-1.
  check(paste0('ordinal-',categories),'nested quadratic',weighted$statistic,ref0$mi*(nrow(d)-1)/nrow(d),2e-4)
  check(paste0('ordinal-',categories),'nested spectrum',weighted$eigenvalues,ref0$scaling.factor*(nrow(d)-1)/nrow(d),2e-4)
}
out<-do.call(rbind,rows)
dir.create('workspace',showWarnings=FALSE)
write.csv(out,'workspace/magmaan-criterion.csv',row.names=FALSE)
write.csv(do.call(rbind,diagnostics),'workspace/magmaan-first-stage-conventions.csv',row.names=FALSE)
files<-list.files(system.file(package='magmaan'),pattern='\\.(rdb|so)$',recursive=TRUE,full.names=TRUE)
capture.output(sessionInfo(),tools::md5sum(files),file='workspace/magmaan-criterion-session.txt')
print(out,row.names=FALSE)
