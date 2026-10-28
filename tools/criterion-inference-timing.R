#!/usr/bin/env Rscript
# Run with OPENBLAS_NUM_THREADS=1 for interpretable small-matrix timings.
pkgload::load_all('.',quiet=TRUE)
model<-'visual =~ x1+x2+x3
textual =~ x4+x5+x6
speed =~ x7+x8+x9'
d<-lavaan::HolzingerSwineford1939
reps<-as.integer(Sys.getenv('INFERENCE_TIMING_REPS','3'))
rows<-list()
for(missing in c(FALSE,TRUE)) {
  set.seed(192)
  x<-d
  if(missing) x$x2[sample(nrow(x),70)]<-NA
  f<-lavaan::cfa(model,x,meanstructure=TRUE,missing=if(missing) 'ml' else 'listwise')
  B<-model_parameter_basis(f)
  theta<-as.numeric(lavaan::coef(f))
  records<-gaussian_pattern_parts(f)$records
  slow<-function(x) gaussian_parameter_evaluation(inference_model_at(f,x),B)$score
  fast<-function(x) gaussian_parameter_gradient(inference_model_at(f,x),B,records)
  ref<-criterion_hessian(slow,theta,B,1e-4)$matrix
  actual<-criterion_hessian(fast,theta,B,1e-4)$matrix
  gap<-max(abs(ref-actual))
  stopifnot(gap<1e-5)
  time<-function(fun) median(replicate(reps,system.time(fun())[['elapsed']]))
  old<-time(function() criterion_hessian(slow,theta,B,1e-4))
  new<-time(function() criterion_hessian(fast,theta,B,1e-4))
  rows[[length(rows)+1L]]<-data.frame(missing,parameters=length(theta),
    casewise_seconds=old,aggregate_seconds=new,speedup=old/new,max_gap=gap)
}
out<-do.call(rbind,rows)
dir.create('workspace',showWarnings=FALSE)
write.csv(out,'workspace/criterion-inference-timing.csv',row.names=FALSE)
capture.output(sessionInfo(),Sys.getenv(c('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')),
               file='workspace/criterion-inference-timing-session.txt')
print(out,row.names=FALSE)
