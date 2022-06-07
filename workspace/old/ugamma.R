ugamma_non_nested <- function(object) {

  # We presently do not support restrictions
  lavmodel <- object@Model
  ceq_idx <- attr(lavmodel@con.jac, "ceq_idx")
  if (length(ceq_idx) > 0L) {
    if (object@SampleStats@ngroups == 1) {
      return(lavInspect(object, "ugamma"))
    } else {
      warning("Restrictions uses 'nested with saturated', not exact match.")
      saturated <- get_saturated(object)
      return(ugamma_nested(object, saturated))
    }
  }

  test <- list()
  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model
  lavoptions <- object@Options
  lavimplied <- object@implied
  lavdata <- object@Data
  test$standard <- object@test[[1]]

  if (test$standard$df == 0L || test$standard$df < 0) {
    stop("Df must be > 0.")
  }

  e <- lavaan:::lav_model_information(
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    extra = TRUE
  )

  delta <- attr(e, "Delta")
  wls_v <- attr(e, "WLS.V")

  gamma <- lavsamplestats@NACOV
  if (is.null(gamma[[1]])) {
    gamma <- lavaan:::lavgamma(object)
  }

  gamma_global <- as.matrix(bdiag(gamma))
  delta_global <- do.call(rbind, delta)
  v_global <- as.matrix(bdiag(wls_v))
  x <- v_global %*% delta_global
  u_global <- v_global - crossprod(t(x), solve(t(delta_global) %*% x, t(x)))
  u_global %*% gamma_global
}

ugamma_nested <- function(m0, m1, method = "delta", a = NULL) {
  # extract information from m1 and m2
  t1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df

  t0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df

  # m = difference between the df's
  m <- r0 - r1

  # check for identical df setting
  if (m == 0L) {
    return(list(
      T.delta = (t0 - t1), scaling.factor = as.numeric(NA),
      df.delta = m, a = as.numeric(NA), b = as.numeric(NA)
    ))
  }

  gamma <- lavTech(m1, "gamma") # the same for m1 and m0
  # check for NULL
  if (is.null(gamma)) {
    stop("lavaan ERROR: Can not compute gamma matrix; perhaps missing \"ml\"?")
  }


  wls_v <- lavTech(m1, "WLS.V")
  pi <- lavaan:::computeDelta(m1@Model)
  p_inv <- lavaan:::lav_model_information_augment_invert(m1@Model,
    information = lavTech(m1, "information"),
    inverted = TRUE
  )

  # compute 'a' matrix
  # NOTE: order of parameters may change between H1 and H0, so be careful!
  if (is.null(a)) {
    a <- lavaan:::lav_test_diff_A(m1, m0, method = method, reference = "H1")
    # take into account equality constraints m1
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  paapaap <- p_inv %*% t(a) %*% MASS::ginv(a %*% p_inv %*% t(a)) %*% a %*% p_inv

  # compute scaling factor
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

  # We need the global gamma, cf. eq.~(10)
  gamma_rescaled <- gamma
  for (i in (seq_along(gamma))) {
    gamma_rescaled[[i]] <- fg[i] * gamma_rescaled[[i]]
  }
  gamma_global <- as.matrix(bdiag(gamma_rescaled))
  # Also the global V:
  v_global <- as.matrix(bdiag(wls_v))
  pi_global <- do.call(rbind, pi)
  # U global version, eq.~(22) in Satorra (2000).
  u_global <- v_global %*% pi_global %*% paapaap %*% t(pi_global) %*% v_global

  return(u_global %*% gamma_global)
}
