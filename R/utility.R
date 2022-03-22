#' Split vector `x` into `n` chunks of equal size.
#' @param x Input vector.
#' @param n Desired output length.
#' @return List of `n` vectors.
#' @keywords internal
chunk <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

#' Calculate distances between an vector and the uniform distribution.
#' @param x Vector of observations in `[0,1]`.
#' @param dist Distance measure.
#' @return The calculated distance.

distance = function(x, dist = c("kolmogorov-smirnov", "anderson-darling")) {
  dist = match.arg(dist)
  if (dist == "kolmogorov-smirnov") {
    n = length(x)
    max(abs(sort(x) - (0:(n - 1)) / n))
  } else {
    stop()
  }
}
