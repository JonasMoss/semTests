function (lavmodel = NULL, lavsamplestats = NULL, lavdata = NULL,
          lavimplied = NULL, lavh1 = NULL, Delta = NULL, lavcache = NULL,
          lavoptions = NULL, extra = FALSE, augmented = FALSE, inverted = FALSE,
          use.ginv = FALSE)
{
  if (.hasSlot(lavmodel, "estimator")) {
    estimator <- lavmodel@estimator
  }
  else {
    estmator <- lavoptions$estimator
  }
  information <- lavoptions$information[1]
  if (is.null(lavh1)) {
    lavh1 <- lav_h1_implied_logl(lavdata = lavdata, lavsamplestats = lavsamplestats,
                                 lavoptions = lavoptions)
  }
  if (information == "observed") {
    if (lavsamplestats@missing.flag || lavdata@nlevels >
        1L) {
      group.weight <- FALSE
    }
    else {
      group.weight <- TRUE
    }
    E <- lav_model_information_observed(lavmodel = lavmodel,
                                        lavsamplestats = lavsamplestats, lavdata = lavdata,
                                        lavimplied = lavimplied, lavh1 = lavh1, lavcache = lavcache,
                                        group.weight = group.weight, lavoptions = lavoptions,
                                        extra = extra, augmented = augmented, inverted = inverted,
                                        use.ginv = use.ginv)
  }
  else if (information == "expected") {
    E <- lav_model_information_expected(lavmodel = lavmodel,
                                        lavsamplestats = lavsamplestats, lavdata = lavdata,
                                        lavimplied = lavimplied, lavh1 = lavh1, lavcache = lavcache,
                                        lavoptions = lavoptions, extra = extra, augmented = augmented,
                                        inverted = inverted, use.ginv = use.ginv)
  }
  else if (information == "first.order") {
    E <- lav_model_information_firstorder(lavmodel = lavmodel,
                                          lavsamplestats = lavsamplestats, lavdata = lavdata,
                                          lavimplied = lavimplied, lavh1 = lavh1, lavcache = lavcache,
                                          lavoptions = lavoptions, check.pd = FALSE, augmented = augmented,
                                          inverted = inverted, use.ginv = use.ginv)
  }
  E
}
