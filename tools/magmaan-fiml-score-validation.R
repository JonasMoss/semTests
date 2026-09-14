#!/usr/bin/env Rscript
stopifnot(requireNamespace("magmaan",quietly=TRUE))
pkgload::load_all(".",quiet=TRUE)
set.seed(907)
d <- lavaan::HolzingerSwineford1939
model <- 'visual =~ x1+x2+x3
textual =~ x4+x5+x6
speed =~ x7+x8+x9'
d$x2[sample(nrow(d),70)] <- NA
d$x5[sample(nrow(d),50)] <- NA
d$x9[d$x1>median(d$x1) & runif(nrow(d))<.4] <- NA
rows <- list()
for(grouped in c(FALSE,TRUE)) {
  args <- list(model=model,data=d,missing="ml",meanstructure=TRUE,fixed.x=FALSE)
  if(grouped) args <- c(args,list(group="school",group.equal="loadings"))
  fit <- do.call(lavaan::cfa,args)
  spec <- magmaan::model_spec(model,meanstructure=TRUE,fixed_x=FALSE,
    group=if(grouped) "school" else NULL,
    group_equal=if(grouped) "loadings" else NULL,
    group_labels=if(grouped) unique(as.character(d$school)) else NULL)
  mag <- magmaan::magmaan(spec,d,estimator="FIML")
  for(sensitivity in c("expected","observed")) for(metric in c("expected","observed")) {
    ref <- tryCatch(magmaan::global_score_flip_test(mag,n_flips=1L,sensitivity=sensitivity,metric=metric),error=identity)
    if(inherits(ref,"error")) {
      out <- tryCatch(score_test(fit,curvature =sensitivity,metric=if(metric=="expected") "information" else "observed"),error=identity)
      stopifnot(grepl("positive definite",conditionMessage(ref)),
                inherits(out,"semTests_error_unstable_spectrum"),
                grepl("positive definite",conditionMessage(out)))
      rows[[length(rows)+1L]] <- data.frame(grouped,sensitivity,metric,
        quantity="both reject indefinite metric",gap=0,tolerance=0)
      next
    }
    out <- score_test(fit,curvature =sensitivity,metric=if(metric=="expected") "information" else "observed")
    full <- score_test(fit,curvature =sensitivity,metric="sandwich")
    values <- list(df=c(out$df,ref$df),statistic=c(out$statistic,ref$statistic_effective),
      spectrum=cbind(sort(out$eigenvalues),sort(ref$eigenvalues)),
      sandwich=c(full$statistic,ref$statistic_sandwich),p=c(out$p_value,ref$p_mixture))
    for(nm in names(values)) {
      v <- matrix(values[[nm]],ncol=2)
      gap <- max(abs(v[,1]-v[,2]))
      tol <- if(nm=="df") 0 else if(nm=="p") 2e-4 else 3e-3
      rows[[length(rows)+1]] <- data.frame(grouped,sensitivity,metric,quantity=nm,gap,tolerance=tol)
      if(gap>tol) stop(paste(grouped,sensitivity,nm,gap))
    }
  }
}
dir.create("workspace",showWarnings=FALSE)
out<-do.call(rbind,rows)
write.csv(out,"workspace/magmaan-fiml-score.csv",row.names=FALSE)
files<-list.files(system.file(package="magmaan"),pattern="\\.(rdb|so)$",recursive=TRUE,full.names=TRUE)
capture.output(sessionInfo(),tools::md5sum(files),file="workspace/magmaan-fiml-score-session.txt")
print(out,row.names=FALSE)
