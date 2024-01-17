
#' @rdname pvalue_internal
pvalues_two <- \(m0, m1, unbiased = FALSE, trad, eba) {
  if (m0@Options$estimator != "ML" || m1@Options$estimator != "ML" ||
      m0@Options$se == "standard" || m1@Options$se == "standard") {
    stop("Only the 'ML' estimator has currently tested.")
  }

  chisq <- lavaan::fitmeasures(m0, "chisq") - lavaan::fitmeasures(m1, "chisq")
  ug <- ugamma_nested(m0, m1, unbiased = unbiased)
  df <- lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]

  peba <- sapply(eba, \(j) eba_pvalue(chisq, lambdas, j))
  names(peba) <- paste0("peba", eba)
  ptrad <- sapply(trad, \(x) trad_pvalue(df, chisq, lambdas, x))
  names(ptrad) <- trad

  c(ptrad, peba)
}


#' Calculate non-nested ugamma for multiple groups.
#' @param object A `lavaan` object.
#' @param unbiased If `TRUE`, uses the unbiased gamma estimate.
#' @keywords internal
#' @return Ugamma for non-nested object.
ugamma_non_nested <- \(object, unbiased = FALSE) {
  lavmodel <- object@Model
  object@Options$gamma.unbiased <- unbiased

  # We do not support restriction fully.
  ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
  if (length(ceq_idx) > 0L) {
    stop("Testing of models with groups and equality constraints not supported.")
    return(ugamma_nested(object, get_saturated(object)))
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

  gamma <- lavaan:::lav_object_gamma(object)
  if (is.null(gamma[[1]])) {
    gamma <- lapply(lavaan::lavInspect(object, "gamma"), \(x) {
      class(x) <- "matrix"
      x
    })
  }

  gamma_global <- as.matrix(Matrix::bdiag(gamma))
  delta_global <- do.call(rbind, delta)
  v_global <- as.matrix(Matrix::bdiag(wls_v))
  x <- v_global %*% delta_global
  u_global <- v_global - crossprod(t(x), solve(t(delta_global) %*% x, t(x)))
  u_global %*% gamma_global
}

#' Calculate nested ugamma.
#'
#' This can also be used with restrictions.
#'
#' @param m0,m1 Two nested `lavaan` objects.
#' @param a The `A` matrix. If if `NULL`, gets calculated by
#'    `lavaan:::lav_test_diff_A` with `method = method`.
#' @param method Method passed to `lavaan:::lav_test_diff_A`.
#' @param unbiased If `TRUE`, uses the unbiased gamma estimate.
#' @keywords internal
#' @return Ugamma for nested object.
ugamma_nested <- \(m0, m1, a = NULL, method = "delta", unbiased = FALSE) {
  m0@Options$gamma.unbiased <- unbiased
  m1@Options$gamma.unbiased <- unbiased

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

  gamma <- lavaan:::lav_object_gamma(m0) # the same for m1 and m0

  # check for NULL
  if (is.null(gamma)) {
    stop("lavaan error: Can not compute gamma matrix; perhaps missing \"ml\"?")
  }

  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")
  # p_inv <- lavaan:::lav_model_information_augment_invert(m1@Model,
  #   information = lavaan::lavTech(m1, "information"),
  #   inverted = TRUE
  # )

  p_inv <- lavaan::lavInspect(m1, what = "inverted.information")

  if (is.null(a)) {
    a <- lavaan:::lav_test_diff_A(m1, m0, method = method, reference = "H1")
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  paapaap <- p_inv %*% t(a) %*% MASS::ginv(a %*% p_inv %*% t(a)) %*% a %*% p_inv

  # compute scaling factor
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

  # We need the global gamma, cf. eq.~(10) in Satorra (2000).
  gamma_rescaled <- gamma
  for (i in (seq_along(gamma))) {
    gamma_rescaled[[i]] <- fg[i] * gamma_rescaled[[i]]
  }
  gamma_global <- as.matrix(Matrix::bdiag(gamma_rescaled))
  # Also the global V:
  v_global <- as.matrix(Matrix::bdiag(wls_v))
  pi_global <- do.call(rbind, pi)
  # U global version, eq.~(22) in Satorra (2000).
  u_global <- v_global %*% pi_global %*% paapaap %*% t(pi_global) %*% v_global
  return(u_global %*% gamma_global)
}
