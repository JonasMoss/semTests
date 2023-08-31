#' Bootstrap `lavaan` objects using the
#'
#' @param ... Several `lavaan` objects.
#' @param functional The functional to calculate. Takes a list of `lavaan`
#'   objects as its argument. Defaults to `identity`, but that's not
#'   a good idea to use, as it consumes a lot of memory.
#' @param n_reps Number of bootstrap repetitions.
#' @param skip_warning If `TRUE`, ignores bootstrapped estimates with
#'   warnings.
#' @keywords internal
#' @return Bootstrapped objects as calculated by `functional`.
bootstrapper <- function(m0, m1 = NULL, functional = identity, n_reps = 1000,
                         iter_max = 100,
                         skip_warning = FALSE) {
  progress <- progressr::progressor(n_reps)

  data <- bollen_stine_transform(m0)
  errors <- 0 # Not in use for the moment.
  future.apply::future_replicate(n_reps,
    {
      result <- NULL
      if (skip_warning) {
        while (is.null(result)) {
          result <- tryCatch(
            {
              boots <- bootstrap(m0, m1, data, iter_max) # lavaan object på bootstrap.
              functional(boots)                # p-values function på lavaan.
            },
            error = function(e) {
              message(paste0("Skipping simulation due to: ", e))
              NULL
            }
          )

          progress()
        }
      } else {
        while (is.null(result)) {
          result <- tryCatch(
            {
              boots <- bootstrap(m0, m1, data)
              functional(boots)
            },
            error = function(e) {
              message(paste0("Skipping simulation due to: ", e))
              NULL
            },
            warning = function(w) {
              message(paste0("Skipping simulation due to: ", w))
              NULL
            }
          )
          progress()
        }
      }
      result
    },
    future.seed = TRUE,
    future.conditions = "message"
  )
}

#' Bootstrap `lavaan` models.
#'
#' @keywords internal
#' @param m0,m1 `lavaan` objects. Data is sample from `m0` and fitted with
#'  `m0` and `m1`. If `m1` is `NULL`, only `m0` is fitted.
#' @param data The data used to sample from, e.g. Bollen-Stine transformed
#'    data.
#' @return A bootstrapped `lavaan` object.
bootstrap <- function(m0, m1 = NULL, data, iter_max = 100) {
  ns <- m0@Data@nobs
  ids <- lapply(ns, function(n) sample(x = n, size = n, replace = TRUE))

  boot_sample <- lavaan::lav_data_update(
    lavdata = m0@Data,
    newX = lapply(seq_along(ns), function(i) data[[i]][ids[[i]], ]),
    lavoptions = lavaan::lavInspect(m0, "options")
  )

  boot_m0 <- lavaan::lavaan(
    slotOptions = m0@Options,
    slotParTable = m0@ParTable,
    slotData = boot_sample,
    control = list(iter.max = iter_max)
  )

  # stopifnot(lavaan::inspect(boot_m0, "converged"))

  if (!is.null(m1)) {
    boot_m1 <- lavaan::lavaan(
      slotOptions = m1@Options,
      slotParTable = m1@ParTable,
      slotData = boot_sample
    )
    stopifnot(lavaan::inspect(boot_m1, "converged"))
    list(boot_m0, boot_m1)
  } else {
    boot_m0
  }
}

#' Bollen-Stine transformer
#'
#' @keywords internal
#' @param object A `lavaan` object.
#' @return A list of Bollen-Stine transformed data.
bollen_stine_transform <- function(object) {
  s <- s_and_s_inv(object)

  lapply(seq(object@SampleStats@ngroups), function(i) {
    data <- object@Data@X[[i]]
    s_sqrt <- s[[i]]$s_sqrt
    s_inv_sqrt <- s[[i]]$s_inv_sqrt
    frame <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% s_sqrt)
    colnames(frame) <- object@Data@ov.names[[i]]
    frame
  })
}

#' Calculate s and s_inv for all subgroups of a `lavaan` object.
#'
#' @param object A `lavaan` object.
#' @keywords internal
#' @return A list containing s and s_inv for all subgroups of a `lavaan` object.
s_and_s_inv <- function(object) {
  lapply(seq(object@SampleStats@ngroups), function(i) {
    s_hat <- lavaan::lav_model_implied(object@Model)$cov[[i]]
    s_inv_hat <- object@SampleStats@icov[[i]]
    list(
      s_sqrt = lavaan::lav_matrix_symmetric_sqrt(s_hat),
      s_inv_sqrt = lavaan::lav_matrix_symmetric_sqrt(s_inv_hat)
    )
  })
}
