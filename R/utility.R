#' Split vector `x` into `n` chunks of equal size.
#' @param x Input vector.
#' @param n Desired output length.
#' @return List of `n` vectors.
#' @keywords internal
chunk <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

#' Estimate the distance between a sample and the uniform distribution.
#'
#' Estimate the distance between a sample and the uniform distribution using
#'   a statistical distance measure.
#'
#' The option `kullback-leibler` uses the Vasicek estimator of the differential
#'   entropy, implemented in `vsgoftest.` The option `cramer-von mises`
#'   returns the Cramer-von Mises distance, the option `anderson-darling`
#'   returns the Anderson-Darling distance (implemented in `goftest`), and
#'   the `kolmogorov-smirnov` option return the Kolmogorov-Smirmov distance. The
#'   option `0.05-distance` measures the distances between the observed
#'   proportion below `0.05` and `0.05` itself.
#'
#' @param x a numeric vector of observations in `[0,1]`.
#' @param dist a distance measure.
#' @export
#' @return Estimated distance between the distribution of `x` and the uniform
#'   distribution.
#'
#' @references
#' Anderson, T.W. and Darling, D.A. (1954) A test of goodness of fit.
#'   Journal of the American Statistical Association 49, 765-69.
#'
#' Cramer, H. (1928). "On the Composition of Elementary Errors".
#'   Scandinavian Actuarial Journal. 1928 (1): 13-4.
#'
#' Vasicek, O., A test for normality based on sample entropy,
#'   Journal of the Royal Statistical Society, 38(1), 54-59 (1976).
#'
#' Daniel, Wayne W. (1990). "Kolmogorov-Smirnov one-sample test".
#'   Applied Nonparametric Statistics (2nd ed.). Boston: PWS-Kent. pp. 319-30.

distance <- function(x, dist = c(
                       "kolmogorov-smirnov",
                       "anderson-darling",
                       "kullback-leibler",
                       "cramer-von mises",
                       "0.05-distance"
                     )) {
  dist <- match.arg(dist)
  stopifnot(all(x >= 0) && all(x <= 1))
  unname(if (dist == "kolmogorov-smirnov") {
    n <- length(x)
    max(abs(sort(x) - (0:(n - 1)) / n))
  } else if (dist == "anderson-darling") {
    goftest::ad.test(x, null = "punif")$statistic / n
  } else if (dist == "cramer-von mises") {
    n <- length(x)
    (1 / (12 * n) + sum(((2 * seq(n) - 1) / (2 * n) - sort(x))^2)) / n
  } else if (dist == "kullback-leibler") {
    vsgoftest::vs.test(x, "dunif", simulate.p.value = FALSE)$statistic
  } else if (dist == "0.05-distance") {
    abs(mean(x <= 0.05) - 0.05)
  })
}

#' Calculate a saturated model.
#' @keywords internal
#' @param object A `lavaan` object.
#' @return A fitted saturated model.
get_saturated <- function(object) {
  data <- lavInspect(object, "data", drop.list.single.group = T)
  data_g <- lapply(seq_len(data), function(x) {
    data_g <- data.frame(data[[x]])
    data_g$g <- x
    data_g
  })
  data <- do.call(rbind, data_g)
  vars <- lavNames(object)
  model <- NULL
  for (i in 1:(length(vars) - 1)) {
    ind <- vars[(i + 1):length(vars)]
    model <- paste(model, ";", paste(vars[i], "~~", paste(ind, collapse = "+")))
  }
  lavaan::cfa(model, data, group = "g", estimator = "MLM")
}

#' Calculate non-nested ugamma.
#' @param object A `lavaan` object.
#' @keywords internal
#' @return Ugamma for non-nested object.
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

#' Calculate non-nested ugamma.
#' @param m0,m1 Two nested `lavaan` objects.
#' @param a The `A` matrix. If if `NULL`, gets calculated by
#'    `lavaan:::lav_test_diff_A` with `method = method`.
#' @param method Method passed to `lavaan:::lav_test_diff_A`.
#' @keywords internal
#' @return Ugamma for non-nested object.
#'
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
