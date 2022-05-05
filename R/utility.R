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
#' @keywords internal
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
  data <- lavaan::lavInspect(object, "data", drop.list.single.group = T)
  data_g <- lapply(seq_along(data), function(x) {
    data_g <- data.frame(data[[x]])
    data_g$g <- x
    data_g
  })
  data <- do.call(rbind, data_g)
  vars <- lavaan::lavNames(object)
  model <- NULL
  for (i in 1:(length(vars) - 1)) {
    ind <- vars[(i + 1):length(vars)]
    model <- paste(model, ";", paste(vars[i], "~~", paste(ind, collapse = "+")))
  }
  estimator <- lavaan::lavInspect(object, "options")$estimator
  lavaan::lavaan(model, data, group = "g", estimator = estimator)
}
