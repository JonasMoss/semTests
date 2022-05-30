#' Bootstrap `lavaan` objects using the
#'
#' @param ... Several `lavaan` objects.
#' @param functional The functional to calculate. Takes a list of `lavaan`
#'   objects as its argument. Defaults to `identity`, but that's not
#'   a good idea to use, as it consumes a lot of memory.
#' @param n_reps Number of bootstrap repetitions.
#' @keywords internal
#' @return Bootstrapped objects as calculated by `functional`.
bootstrapper <- function(..., functional = identity, n_reps = 1000) {
  progress <- progressr::progressor(n_reps)

  data <- bollen_stine_transform(...)
  errors <- 0 # Not in use for the moment.
  models <- list(...)
  future.apply::future_replicate(n_reps, {
    result <- NULL
    while (is.null(result)) {
      result <- tryCatch({
          boots <- lapply(seq_along(models), function(i) {
            bootstrap(models[[i]], data[[i]])
          })
          functional(boots)
        },
        error = function(e) {
          print(e)
          errors <<- errors + 1
          NULL
        }
      )
      progress()
    }
    result
  }, future.seed = TRUE)
}

#' Bootstrap a single model with groups one time.
#'
#' @keywords internal
#' @param object A `lavaan` object.
#' @param data The data used to sample from, e.g. Bollen-Stein transformed
#'    data.
#' @return A bootstrapped `lavaan` object.
bootstrap <- function(object, data) {
  ns <- object@Data@nobs
  ids <- lapply(ns, function(n) sample(x = n, size = n, replace = TRUE))

  boot_sample <- lavaan::lav_data_update(
    lavdata = object@Data,
    newX = lapply(seq_along(ns), function(i) data[[i]][ids[[i]], ]),
    lavoptions = lavaan::lavInspect(object, "options")
  )

  boot_object <- lavaan::lavaan(
    slotOptions = object@Options,
    slotParTable = object@ParTable,
    slotData = boot_sample
  )

  stopifnot(lavaan::inspect(boot_object, "converged"))
  boot_object
}

#' Bollen-Stine transformer
#'
#' Does the Bollen-Stine transform on a list of models using the implied
#'    covariance structure of the first model.
#'
#' @keywords internal
#' @param ... Models to transform.
#' @return A list of Bollen-Stine transformed data.
bollen_stine_transform <- function(...) {
  models <- list(...)

  # Calculate s_hat and s_inv_hat for all groups of the first model.
  s <- s_and_s_inv(models[[1]])

  # Transform the data for each model.
  lapply(models, function(model) {
    lapply(seq(model@SampleStats@ngroups), function(i) {
      data <- model@Data@X[[i]]
      s_sqrt <- s[[i]]$s_sqrt
      s_inv_sqrt <- s[[i]]$s_inv_sqrt
      frame <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% s_sqrt)
      colnames(frame) <- model@Data@ov.names[[i]]
      frame
    })
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
