hs_model <- " visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 "

## Estimation that IS allowed.
object <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939, estimator = "MLM", test = "sb"
)

m1 <- lavaan::cfa(hs_model,
                  data = lavaan::HolzingerSwineford1939,
                  group = "school", estimator = "ML", test = "satorra.bentler"
)

unbiased <- FALSE
collapse <- FALSE

stopifnot(is.logical(unbiased))
lavoptions = lavaan::lavInspect(object, "options")
lavdata = object@Data
lavoptions$gamma.unbiased <- unbiased
gamma_list <- lavaan:::lav_samplestats_from_data(lavdata, lavoptions)@NACOV

pvalues(object, tests = "SB_UG")



