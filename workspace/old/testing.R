lavaan:::lav_test_satorra_bentler(object, test = "scaled.shifted")



df <- object@test$standard$df
ug <- lavaan::inspect(object, "UG")
trace_ug <- sum(diag(ug))
trace_ug2 <- sum(diag(ug %*% ug))
group <- object@test$standard$stat.group
lavsamplestats = object@SampleStats
fg <- unlist(lavsamplestats@nobs)/lavsamplestats@ntotal
a <- sqrt(df/trace_ug2)
shift.parameter <- fg * (df - a*trace_ug)
scaling_factor  <- 1/a
if(scaling_factor < 0) scaling_factor <- as.numeric(NA)
stat_group <- (group * a + shift.parameter)
stat <- sum(stat_group)
1 - pchisq(stat, df)

### mean-var adjusted
trace_ug <- sum(diag(ug))
trace_ug2 <- sum(diag(ug %*% ug))
df <- trace_ug ^ 2 / trace_ug2
scaling_factor <- trace_ug/df
group <- object@test$standard$stat.group
if(scaling_factor < 0) scaling_factor <- as.numeric(NA)
stat_group <- group/ scaling_factor
stat <- sum(stat_group)
1 - pchisq(stat, df)


lavaan:::lav_test_satorra_bentler(object, test = "mean.var.adjusted")
