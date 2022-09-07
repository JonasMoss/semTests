#' Select *p*-value for a lavaan model using the bootstrap selector.
#'
#' @param m0,m1 One or two `lavaan` objects. If two, the first object should be
#'    with restrictions and the second without.
#' @param n_reps Number of bootstrap repetitions.
#' @param distances A vector of strings containing the distances to calculate.
#'    Passed to `distance`.
#' @param skip_warning If `TRUE`, ignores bootstrapped estimates with
#'   warnings.
#' @export
#' @return An object of class `semselector`, inheriting from `data.frame`,
#'    containing the winning *p*-values for the selected distances.
#'    The attributes are:
#'    * `boots`: The data frame of bootstrap samples.
#'    * `n_reps`: Argument passed to `n_reps`.
#'    * `pvalues`: The calculated p-values.
#'    * `distances` The estimated distances.
#'    * `bollen-stine`: The Bollen-Stine *p*-value.

semselector <- function(m0, m1 = NULL,
                        n_reps = 1000,
                        distances = c(
                          "kolmogorov-smirnov",
                          "anderson-darling",
                          "cramer-von mises",
                          "kullback-leibler",
                          "0.05-distance"
                        ),
                        skip_warning = FALSE) {
  pvals <- pvalues(m0, m1)

  samples <- if (!is.null(m1)) {
    bootstrapper(
      m0,
      m1,
      functional = function(x) pvalues_two(x[[1]], x[[2]]),
      n_reps = n_reps,
      skip_warning = skip_warning
    )
  } else {
    bootstrapper(
      m0,
      functional = function(x) pvalues_one(x),
      n_reps = n_reps,
      skip_warning = skip_warning
    )
  }


  samples <- pmin(pmax(samples, 0), 1)

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
  bs <- mean(attr(minimals, "boots")$pstd < pvals[1])
  attr(minimals, "bollen-stine") <- bs
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
