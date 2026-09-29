#!/usr/bin/env Rscript
# Reproducible calibration preflight, not a paper simulation bank.
# All attempts, including fitting/curvature failures, are retained.
pkgload::load_all(".",quiet=TRUE)
reps<-as.integer(Sys.getenv("INFERENCE_REPS","200"))
stopifnot(is.finite(reps),reps>0)
set.seed(420190)
rows<-list()
record<-function(design,rep,fun) {
  result<-tryCatch(suppressWarnings(fun()),error=function(e) list(error=conditionMessage(e)))
  if(is.list(result) && !is.null(result$error)) {
    return(data.frame(design,rep,method="all",p=NA_real_,error=result$error))
  }
  data.frame(design,rep,method=names(result),p=as.numeric(result),error="")
}
quality<-function(...) for(f in list(...)) {
  if(!isTRUE(lavaan::lavInspect(f,"converged")) || !isTRUE(lavaan::lavInspect(f,"post.check")))
    stop("Fit failed convergence or admissibility check")
}
# E[X]=1, E[Y]=.5, Var(X)=Var(Y)=1, Cov(X,Y)=.3.
# Both fitted means are fixed at zero: beta*=E[XY]/E[X^2]=.8/2=.4.
for(i in seq_len(reps)) {
  rows[[length(rows)+1L]]<-record("ML misspecified means, beta*=0.4",i,function() {
    z<-matrix(rnorm(1200),600,2)
    d<-data.frame(x=1+z[,1],y=.5+.3*z[,1]+sqrt(.91)*z[,2])
    model<-'y ~ b*x\nx ~ 0*1\ny ~ 0*1'
    f<-lavaan::sem(model,d,fixed.x=FALSE,meanstructure=TRUE)
    f0<-lavaan::sem(paste(model,'b==.4',sep='\n'),d,fixed.x=FALSE,meanstructure=TRUE)
    quality(f,f0)
    s<-score_test_nested(f0,f,keep=TRUE)
    full<-score_test_matrices(s$ingredients$score,NULL,s$ingredients$score_covariance)
    stat<-2*s$nobs*(lavaan::fitMeasures(f0,'fmin')-lavaan::fitMeasures(f,'fmin'))
    c(score=full$p_value,curvature_score=s$p_value,
      criterion=pvalue_all(as.numeric(stat),s$eigenvalues))
  })
}
# Two independent identical ordinal populations imply equal pseudo-loadings
# under the same DWLS criterion, despite a misspecified one-factor structure.
for(i in seq_len(reps)) {
  rows[[length(rows)+1L]]<-record("DWLS common pseudo-loadings, two identical groups",i,function() {
    n<-800
    z<-matrix(rnorm(n*6),n,6)
    f1<-z[,1];f2<-.5*f1+sqrt(.75)*z[,2]
    x<-cbind(.7*f1+sqrt(.51)*z[,3],.6*f1+.8*z[,4],
             .7*f2+sqrt(.51)*z[,5],.6*f2+.8*z[,6])
    d<-as.data.frame(1L*(x>0));names(d)<-paste0('x',1:4)
    d$g<-rep(c('a','b'),each=n/2)
    model<-'f =~ NA*x1+x2+x3+x4\nf ~~ 1*f'
    args<-list(model=model,data=d,ordered=paste0('x',1:4),group='g')
    f<-do.call(lavaan::cfa,args)
    f0<-do.call(lavaan::cfa,c(args,list(group.equal='loadings')))
    quality(f,f0)
    s<-score_test_nested(f0,f,keep=TRUE)
    full<-score_test_matrices(s$ingredients$score,NULL,s$ingredients$score_covariance)
    stat<-2*s$nobs*(lavaan::fitMeasures(f0,'fmin')-lavaan::fitMeasures(f,'fmin'))
    c(score=full$p_value,curvature_score=s$p_value,
      criterion=pvalue_all(as.numeric(stat),s$eigenvalues))
  })
  if(i%%25==0) cat("Ordinal attempts:",i,"/",reps,"\n")
}
out<-do.call(rbind,rows)
dir.create('workspace',showWarnings=FALSE)
write.csv(out,'workspace/inference-calibration-attempts.csv',row.names=FALSE)
summaries<-lapply(split(out[out$error=='',],interaction(out[out$error=='','design'],out[out$error=='','method'],drop=TRUE)),function(x) {
  rate<-mean(x$p<.05)
  data.frame(design=x$design[1],method=x$method[1],attempts=reps,successful=nrow(x),
             rejected=sum(x$p<.05),rejection_rate=rate,mcse=sqrt(rate*(1-rate)/nrow(x)))
})
summary<-do.call(rbind,summaries)
write.csv(summary,'workspace/inference-calibration-summary.csv',row.names=FALSE)
capture.output(sessionInfo(),file='workspace/inference-calibration-session.txt')
print(summary,row.names=FALSE)
cat('Failed attempts:',sum(out$error!=''),'of',2*reps,'; all retained in attempts CSV.\n')
