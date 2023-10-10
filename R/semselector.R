#' Select *p*-value for a lavaan model using the bootstrap selector.
#'
#' @param m0,m1 One or two `lavaan` objects. If two, the first object should be
#'    with restrictions and the second without.
#' @param n_reps Number of bootstrap repetitions.
#' @param distances A vector of strings containing the distances to calculate.
#'    Passed to `distance`.
#' @param trad,eba The choice of p-values to calculate.
#' @param unbiased Use unbiased gamma estimate? 1: No, 2: Yes, 3: Use both.
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
                        trad = list("pstd", "psb", "pss"),
                        eba = c(2, 4),
                        unbiased = 1,
                        chisq = c("trad", "rls"),
                        bollen_stine = c(1, 2, 3),
                        distances = c(
                          "kolmogorov-smirnov",
                          "anderson-darling",
                          "cramer-von mises",
                          "0.05-distance"
                        ),
                        bs = TRUE,
                        do_bs = TRUE,
                        skip_warning = FALSE) {
  pvals <- pvalues(m0, m1, trad = trad, eba = eba, unbiased = unbiased, chisq = chisq)

  samples <- if (!is.null(m1)) {
    bootstrapper(
      m0,
      m1,
      functional = \(x) pvalues(x[[1]], x[[2]], trad = trad, eba = eba, unbiased = unbiased, chisq = chisq),
      n_reps = n_reps,
      skip_warning = skip_warning
    )
  } else {
    bootstrapper(
      m0,
      functional = \(x) pvalues(x, trad = trad, eba = eba, unbiased = unbiased, chisq = chisq),
      n_reps = n_reps,
      bs = bs,
      skip_warning = skip_warning
    )
  }

  samples <- pmin(pmax(samples, 0), 1)

  boot_dists <- sapply(distances, function(d) apply(samples, 1, distance, d))
  output <- data.frame(
    apply(boot_dists, 2, min),
    rownames(boot_dists)[apply(boot_dists, 2, which.min)],
    pvals[apply(boot_dists, 2, which.min)]
  )
  colnames(output) <- c("distance", "type", "pvalue")

  class(output) <- c("semselector", "data.frame")
  attr(output, "boots") <- as.data.frame(t(samples))
  attr(output, "n_reps") <- n_reps
  attr(output, "distances") <- boot_dists

  if ("trad" %in% chisq) {
    if (bs) {
      bs_trad <- mean(attr(output, "boots")$pstd_trad < pvals["pstd_trad"])
    } else {
      f <- \(object) {
        chisq <- lavaan::fitmeasures(object, "chisq")
        df <- lavaan::fitmeasures(object, "df")
        1 - stats::pchisq(chisq, df)
      }
      bss <- bootstrapper(
        m0 = m0, m1 = NULL, functional = f, n_reps = n_reps, bs = TRUE,
        skip_warning = skip_warning
      )
      bs_trad <- mean(bss < pvals["pstd_trad"])
    }
    pvals <- c(pvals, pbs_trad = bs_trad)
  }

  if ("rls" %in% chisq) {
    if (bs) {
      bs_rls <- mean(attr(output, "boots")$pstd_rls < pvals["pstd_rls"])
    } else {
      f <- \(object) {
        chisq <- rls(object)
        df <- lavaan::fitmeasures(object, "df")
        1 - stats::pchisq(chisq, df)
      }
      bss <- bootstrapper(
        m0 = m0, m1 = NULL, functional = f, n_reps = n_reps, bs = TRUE,
        skip_warning = skip_warning
      )
      bs_rls <- mean(bss < pvals["pstd_rls"])
    }
    pvals <- c(pvals, pbs_rls = bs_rls)
  }

  attr(output, "pvalues") <- pvals
  output
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
