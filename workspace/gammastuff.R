obj@Options$gamma.unbiased <- FALSE
gamma <- lavaan:::lav_object_gamma(obj)[[1]]


gamma_unbiased <- gamma_est_unbiased(n = lavaan::lavInspect(obj, "nobs"),
                   sigma = obj@SampleStats@cov[[1]],
                   gamma_adf = gamma,
                   gamma_nt = NULL)

sigma <- obj@SampleStats@cov[[1]]
sum(abs(2 * lavaan:::lav_matrix_duplication_ginv_pre_post(sigma %x% sigma) - gamma_est_nt(sigma)))


sum(abs(tcrossprod(vech(sigma)) - tcrossprod(lavaan:::lav_matrix_vech(sigma))))

obj@Options$gamma.unbiased <- TRUE
gamma_unbiased2 <- lavaan:::lav_object_gamma(obj)[[1]]

sum(abs(gamma_unbiased - gamma_unbiased2))
sum(abs(gamma  - gamma_unbiased2))


obj@SampleStats@cov[[1]] - cov(x[, 1:25], use = "na.or.complete")


n <- nrow(na.omit(x))
lavaan::lavInspect(obj, "sampstat")$cov * n / (n-1) - cov(x[, 1:25], use = "complete.obs")
obj@SampleStats@cov - cov(x[, 1:25], use = "complete.obs")

