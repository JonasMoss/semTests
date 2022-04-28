#' Select p-value for a lavaan model using the bootstrap selector.
#'
#' @param object A fitted `lavaan` object.
#' @param n_reps Number of bootstrap repetitions.
#' @param distances A vector of strings containing the distances to calculate.
#'    Passed to `distance`. **Note:** `kullback-leibler` is very slow, hence not
#'    in the default arguments.
#' @param custom List of named custom functions passed that take
#'   `lavaan` objects as input.
#' @param type Normal parametric bootstrap / transformed samples bootstrap.
#' @export
#' @return An object of class `semselector`.

semselector <- function(object,
                        n_reps = 1000,
                        distances = c(
                          "kolmogorov-smirnov",
                          "anderson-darling",
                          "cramer-von mises",
                          "kullback-leibler",
                          "0.05-distance"
                        ),
                        custom = NULL,
                        type = NULL) {
  pvals <- pvalues(object)
  samples <- bootstrapper(object, n_reps)
  boot_dists <- sapply(distances, function(d) apply(samples, 1, distance, d))
  minimals <- data.frame(
    apply(boot_dists, 2, min),
    rownames(boot_dists)[apply(boot_dists, 2, which.min)],
    pvals[apply(boot_dists, 2, which.min)]
  )
  colnames(minimals) <- c("distance", "type", "pvalue")

  class(minimals) <- c("semselector", "data.frame")
  attr(minimals, "boots") <- as.data.frame(t(samples))
  attr(minimals, "n_reps") <- n_reps
  attr(minimals, "pvalues") <- pvals
  attr(minimals, "distances") <- boot_dists
  minimals
}

#' Plotting generic for `semselector` objects.
#' @param x `semselector` object.
#' @param y Ignored.
#' @param nrow Passed to `ggplot2::facet_wrap`.
#' @param binwidth Passed to `ggplot2::geom_histogram`. Defaults to `0.05`, to
#'   make it easier to see the significance cutoff.
#' @param ... Passed to `ggplot2::geom_histogram`.
#' @export
plot.semselector <- function(x, y, nrow = 3, binwidth = 0.05, ...) {
  value <- NULL # To avoid CRAN check note.
  data <- dplyr::arrange(tidyr::pivot_longer(
    attr(x, "boots"),
    tidyr::everything()
  ), name)
  ggplot2::ggplot(data, ggplot2::aes(value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
      binwidth = binwidth,
      boundary = 0, ...
    ) +
    ggplot2::xlab("p-value") +
    ggplot2::facet_wrap(~name, nrow = nrow)
}
