U <- lavaan::lavInspect(object, "U")
gamma <- lavaan::lavInspect(object, "gamma")
Ugamma_ <- lavaan::lavInspect(object, "Ugamma")

mean(Ugamma_ == U %*% gamma)


library("lavaan")
lavTech(m1, "gamma")



str(lavaan:::lavInspect.lavaan(m1, what = "gamma", add.labels = FALSE,
                  add.class = FALSE, list.by.group = FALSE,
                  drop.list.single.group = FALSE))


str(lavaan:::lav_object_inspect_sampstat_gamma(m1,
                                  add.labels = FALSE,
                                  add.class = FALSE,
                                  drop.list.single.group = FALSE))

lavaan:::lav_object_inspect_sampstat_gamma(m1,
                                           add.labels = FALSE,
                                           add.class = FALSE,
                                           drop.list.single.group = FALSE)


lavaan:::lav_object_gamma(m1)


x <- lavaan::lavTech(m1, "data")


a <- lavaan:::lav_samplestats_Gamma(x[[1]],
                                  Mu                 = NULL,
                                  Sigma              = NULL,
                                  x.idx              = integer(0L),
                                  cluster.idx        = NULL,
                                  fixed.x            = FALSE,
                                  conditional.x      = FALSE,
                                  meanstructure      = TRUE,
                                  slopestructure     = FALSE,
                                  gamma.n.minus.one  = FALSE,
                                  unbiased           = TRUE,
                                  Mplus.WLS          = FALSE)

b <- gamma_est_unbiased(x[[1]])


b <- lavTech(m1, "gamma")[[1]]


caller <- \(x) lavaan:::lav_samplestats_Gamma(x, meanstructure = TRUE, unbiased = TRUE)
meanstructure <- m1@Options$meanstructure
groups <- m1@Data@ngroups
out <- vector("list", length = groups)
for (g in seq_len(groups)) {
  out[[g]] = caller(x[[g]])
}

m1

m1@Options$gamma.unbiased = FALSE
a <- lavaan:::lav_object_gamma(m1)

m1@Options$gamma.unbiased = TRUE
b <- lavaan:::lav_object_gamma(m1)
