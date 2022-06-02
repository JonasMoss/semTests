#' Calculate p-values for one or two lavaan objects.
#'
#' Calculate p-values for a `lavaan` object using several methods.
#'
#' * `pstd` the standard *p*-value extracted from lavaan.
#' * `psb` Satorra-Bentler *p*-value.
#' * `pfull` *p*-value based on all eigenvalues of the gamma matrix.
#' * `phalf` *p*-value based on largest half of eigenvalues of the gamma matrix.
#' * `pcf` Scaled F *p*-value.
#' * `pss` Scaled and shifted *p*-value.
#'
#' @param m0,m1 One or two `lavaan` objects. If two, the first object should be
#'    with restrictions and the second without.
#' @export
#' @return A named vector containing the p-values `pml`. `psb`, `pfull`,
#'    `phalf`, `pcf`, `pss`.
pvalues <- function(m0, m1) {
  if (missing(m1)) {
    pvalues_one(m0)
  } else {
    pvalues_two(m0, m1)
  }
}


#' P value function for one and two arguments.
#'
#' @keywords internal
#' @param object,m0,m1 `lavaan` object.
#' @name pvalue_internal
#' @return pvalues.
NULL

#' @rdname pvalue_internal
pvalues_one <- function(object) {
  if (object@Options$estimator != "ML") {
    stop("Only the 'ML' estimator has currently tested.")
  }

  chisq <- lavaan::fitmeasures(object, "chisq")
  ug <- ugamma_non_nested(object)
  df <- lavaan::fitmeasures(object, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]
  eigenps <- eigen_pvalues(chisq, lambdas)

  c(
    pstd = unname(lavaan::fitmeasures(object, "pvalue")),
    psb = eigenps$psb,
    pfull = eigenps$pfull,
    phalf = eigenps$phalf,
    plog = eigenps$plog,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(object),
    pmv = mean_var_adjusted(object)
  )
}

#' @rdname pvalue_internal
pvalues_two <- function(m0, m1) {
  if (m0@Options$estimator != "ML" || m1@Options$estimator != "ML") {
    stop("Only the 'ML' estimator has currently tested.")
  }

  aov <- lavaan::anova(m1, m0)
  chisq <- aov$`Chisq diff`[[2]]

  ug <- ugamma_nested(m0, m1)
  df <- lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")
  lambdas <- Re(eigen(ug)$values)[seq(df)]
  eigenps <- eigen_pvalues(chisq, lambdas)

  c(
    pstd = aov$`Pr(>Chisq)`[[2]],
    psb = eigenps$psb,
    pfull = eigenps$pfull,
    phalf = eigenps$phalf,
    plog = eigenps$plog,
    psf = scaled_f(chisq, lambdas),
    pss = scaled_and_shifted(m0, m1),
    pmv = mean_var_adjusted(m0, m1)
  )
}

#' Calculate the scaled and shifted / the mean-variance adjusted p-value
#'
#' @param m0,m1 `lavaan` objects.
#' @name laavan_tests
#' @return The scaled and shifted p-value or the mean-variance adjusted p-value.
NULL

#' @rdname laavan_tests
scaled_and_shifted <- function(m0, m1 = NULL) {

  if (is.null(m1)) {
    m <- suppressWarnings(lavaan::lavaan(
      slotOptions = m0@Options,
      slotParTable = m0@ParTable,
      slotData = m0@Data,
      test = "scaled.shifted"
    ))
    unname(lavaan::fitmeasures(m, fit.measures = "pvalue.scaled"))
  } else {
    lavaan::lavTestLRT(m0, m1, method = "satorra.2000", scaled.shifted = TRUE)$`Pr(>Chisq)`[2]
  }
}

#' @rdname laavan_tests
mean_var_adjusted <- function(m0, m1 = NULL) {

  if (is.null(m1)) {
    m <- suppressWarnings(lavaan::lavaan(
      slotOptions = m0@Options,
      slotParTable = m0@ParTable,
      slotData = m0@Data,
      test = "mean.var.adjusted"
    ))
    unname(lavaan::fitmeasures(m, fit.measures = "pvalue.scaled"))
  } else {
    lavaan::lavTestLRT(m0, m1, method = "satorra.2000", scaled.shifted = FALSE)$`Pr(>Chisq)`[2]
  }

}

#' Calculate the scaled_f p-value.
#' @param chisq Chi-square fit value from a lavaan object.
#' @param eig eig of UG matrix.
#' @return scaled f p-value.
#' @keywords internal
scaled_f <- function(chisq, eig) {
  s1 <- sum(eig)
  s2 <- sum(eig^2)
  s3 <- sum(eig^3)

  denom <- 2 * s1 * s2^2 - s1^2 * s3 + 2 * s2 * s3
  if (denom > 0) {
    d1f3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) / denom
    d2f3 <- (s1^2 * s2 + 2 * s2^2) / (s3 * s1 - s2^2) + 6
    if (d2f3 < 6) d2f3 <- Inf
    cf3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) /
      (s1^2 * s2 - 4 * s2^2 + 6 * s1 * s3)
  } else {
    d1f3 <- Inf
    d2f3 <- s1^2 / s2 + 4
    cf3 <- s1 * (s1^2 + 2 * s2) / (s1^2 + 4 * s2)
  }

  unname(1 - stats::pf(chisq / cf3, d1f3, d2f3))
}

#' Get eigenvalue-based p-values.
#' @param chisq Chi-square fit value from a lavaan object.
#' @param eig Eigenvalues of the UG matrix.
#' @return List of eigenvalue-based p-values. `psb` is Satorra--Bentler,
#'    `pfull` is based on every p-value, while `phalf` is based on the
#'    half of the (largest) p-values.
#' @keywords internal
eigen_pvalues <- function(chisq, eig) {
  m <- length(eig)
  pfull <- CompQuadForm::imhof(chisq, eig)$Qq
  psb <- CompQuadForm::imhof(chisq, rep(mean(eig), m))$Qq
  alpha <- 0.0
  lambdas_log <- mean(eig) + (1 - alpha) * (mean(log(seq(m))) - log(seq(m)))
  plog <- CompQuadForm::imhof(chisq, lambdas_log)$Qq

  k <- ceiling(m / 2)
  eig[1:k] <- mean(eig[1:k])
  eig[(k + 1):m] <- mean(eig[(k + 1):m])

  phalf <- CompQuadForm::imhof(chisq, eig)$Qq
  list(pfull = pfull, phalf = phalf, psb = psb, plog = plog)
}


#' Calculate non-nested ugamma for multiple groups.
#' @param object A `lavaan` object.
#' @keywords internal
#' @return Ugamma for non-nested object.
ugamma_non_nested <- function(object) {

  # We presently do not support restrictions
  lavmodel <- object@Model

  if (object@SampleStats@ngroups == 1) {
    return(lavaan::lavInspect(object, "Ugamma"))
  }

  # We presently do not support restriction fully.
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

  gamma <- lavsamplestats@NACOV
  if (is.null(gamma[[1]])) {
    gamma <- lapply(lavaan::lavInspect(object, "gamma"), function(x) {
      class(x) = "matrix"
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
#' @keywords internal
#' @return Ugamma for non-nested object.
ugamma_nested <- function(m0, m1, a = NULL, method = "delta") {
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

  gamma <- lavaan::lavTech(m1, "gamma") # the same for m1 and m0
  # check for NULL
  if (is.null(gamma)) {
    stop("lavaan error: Can not compute gamma matrix; perhaps missing \"ml\"?")
  }


  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")
  p_inv <- lavaan:::lav_model_information_augment_invert(m1@Model,
    information = lavaan::lavTech(m1, "information"),
    inverted = TRUE
  )

  # compute A matrix
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
  gamma_global <- as.matrix(Matrix::bdiag(gamma_rescaled))
  # Also the global V:
  v_global <- as.matrix(Matrix::bdiag(wls_v))
  pi_global <- do.call(rbind, pi)
  # U global version, eq.~(22) in Satorra (2000).
  u_global <- v_global %*% pi_global %*% paapaap %*% t(pi_global) %*% v_global

  return(u_global %*% gamma_global)
}
