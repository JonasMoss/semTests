pvalues(m0, m1, ) {
  if (missing(m1)) {
    pvalues_one_model
  } else {
    pvalues_two_models
  }
}


pvalues_two_model(m0, m1, ) {
  if (missing(m1)) {
    pvalues_one_model
  } else {
    pvalues_two_models
  }
}



f = function(...) {
  print(list(...))
}

g = function(m0, m1) {
  f(m0, m1)
}



f = function(..., n_reps = 1) {
  ..1
}


# generate syntax for an independence model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)

# extract data slot and options
lavdata <- fit@Data
lavoptions <- lavInspect(fit, "options")

# generate sample statistics object
sampleStats <- lav_samplestats_from_data(lavdata = lavdata,
                                         lavoptions = lavoptions)


fit.h0 <- lavaan(slotOptions     = fit@Options,
                 slotParTable    = fit@ParTable,
                 slotData        = fit@Data)
