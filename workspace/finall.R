require(lavaan)
require(Matrix)

UGamma_nonNested <- function(lavobject) {

  #We presently do not support restrictions
  lavmodel <- lavobject@Model
  ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
  if(length(ceq.idx) > 0L) {

    if(lavobject@SampleStats@ngroups == 1) {
      return(lavInspect(lavobject, "UGamma"))
    } else{
      warning("UGamma_nonNested with restrictions uses nested with saturated, not exact match.")
      saturated <- getsaturated(lavobject)
      return(UGamma_nested(lavobject, saturated))
    }
  }

  TEST <- list()
  lavsamplestats <- lavobject@SampleStats
  lavmodel       <- lavobject@Model
  lavoptions     <- lavobject@Options
  lavpartable    <- lavobject@ParTable
  lavimplied     <- lavobject@implied
  lavdata        <- lavobject@Data
  TEST$standard  <- lavobject@test[[1]]

  if (TEST$standard$df == 0L || TEST$standard$df < 0) {
    stop("Df must be > 0.")
  }

  E <- lavaan:::lav_model_information(
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    extra = TRUE
  )
  Delta <- attr(E, "Delta")
  WLS.V <- attr(E, "WLS.V")

  Gamma <- lavsamplestats@NACOV
  if (is.null(Gamma[[1]])) {
    Gamma <- lavaan:::lavGamma(lavobject)
  }

  GammaGlobal <- as.matrix(bdiag(Gamma))
  DeltaGlobal <- do.call(rbind, Delta)
  VGlobal <- as.matrix(bdiag(WLS.V))
  J.inv <- solve(t(DeltaGlobal) %*% VGlobal %*% DeltaGlobal)
  UGlobal <-
    VGlobal - VGlobal %*% DeltaGlobal %*% J.inv %*% t(DeltaGlobal) %*% VGlobal
  #UGammaGlobal <- UGlobal %*% GammaGlobal
  #sum(diag(UGammaGlobal))

  #return(list(UGamma=UGlobal %*% GammaGlobal, df=lavobject))
  UGlobal %*% GammaGlobal
}

UGamma_nested <- function(m0, m1, A.method = "delta", A = NULL){#m0 restricted model.
  # extract information from m1 and m2
  T1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df

  T0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df

  # m = difference between the df's
  m <- r0 - r1

  # check for identical df setting
  if(m == 0L) {
    return(list(T.delta = (T0 - T1), scaling.factor = as.numeric(NA),
                df.delta = m, a = as.numeric(NA), b = as.numeric(NA)))
  }

  Gamma <- lavTech(m1, "Gamma") # the same for m1 and m0
  # check for NULL
  if(is.null(Gamma)) {
    stop("lavaan ERROR: can not compute Gamma matrix; perhaps missing = \"ml\"?")
  }


  WLS.V <- lavTech(m1, "WLS.V")
  PI <- lavaan:::computeDelta(m1@Model)
  P <- lavTech(m1, "information")
  # needed? (yes, if H1 already has eq constraints)
  P.inv <- lavaan:::lav_model_information_augment_invert(m1@Model,
                                                         information = P,
                                                         inverted = TRUE)
  # compute 'A' matrix
  # NOTE: order of parameters may change between H1 and H0, so be
  # careful!
  if(is.null(A)) {
    A <- lavaan:::lav_test_diff_A(m1, m0, method = A.method, reference = "H1")
    # take into account equality constraints m1
    if(m1@Model@eq.constraints) {
      A <- A %*% t(m1@Model@eq.constraints.K)
    }
  }

  # compute tr UG per group
  ngroups <- m1@SampleStats@ngroups
  UG.group  <- vector("list", length=ngroups)

  # safety check: A %*% P.inv %*% t(A) should NOT contain all-zero
  # rows/columns
  # FIXME: is this really needed? As we use ginv later on
  APA <- A %*% P.inv %*% t(A)
  cSums <- colSums(APA)
  rSums <- rowSums(APA)
  empty.idx <- which( abs(cSums) < .Machine$double.eps^0.5 &
                        abs(rSums) < .Machine$double.eps^0.5 )
  if(length(empty.idx) > 0) {
    A <- A[-empty.idx,, drop = FALSE]
  }

  # PAAPAAP
  PAAPAAP <- P.inv %*% t(A) %*% MASS::ginv(A %*% P.inv %*% t(A)) %*% A %*% P.inv

  # compute scaling factor
  fg <- unlist(m1@SampleStats@nobs)/m1@SampleStats@ntotal

  #We need the global gamma, cf. eq.~(10)
  GammaRescaled <- Gamma
  for(i in (1:length(Gamma))) {
    GammaRescaled[[i]] <- fg[i]*GammaRescaled[[i]]
  }
  GammaGlobal <- as.matrix(bdiag(GammaRescaled))
  #Also the global V:
  VGlobal <- as.matrix(bdiag(WLS.V))
  PiGlobal <- do.call(rbind, PI)
  #U global version, eq.~(22) in Satorra (2000).
  UGlobal <- VGlobal %*% PiGlobal %*% PAAPAAP %*% t(PiGlobal) %*% VGlobal

  return(UGlobal %*% GammaGlobal)
}

#brute force. This is only called for two-or-more groups

getsaturated <- function(object){
  # dat = lavInspect(object, "data", drop.list.single.group = T)
  dat = data.frame(as.matrix(object@Data@X))
  names(dat) <- object@Data@ov.names[[1]]
  # ngroups <- length(dat)
  # datg <- lapply(1:ngroups, function(x){
  #   dd <- data.frame(dat[[x]])
  #   dd$g  <- x
  #   dd
  # })
  dat = do.call(rbind, datg)
  vars <- lavNames(object)
  model <- NULL
  for(i in 1:(length(vars)-1)){
    model <- paste(model, ";", paste(vars[i], "~~", paste(vars[(i+1):length(vars)], collapse="+")))
  }

  lavaan::cfa(model, dat, group="g", estimator="MLM")
}
