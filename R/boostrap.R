#' Bootstrap p-values for a `lavaan` object.
#' @param object A `lavaan` object.
#' @param n_reps Number of bootstrap repetitions.
#' @export
#' @return Bootstrapped p-values.

bootstrapper = function(object, n_reps = 10) {

  progress <- progressr::progressor(n_reps)
  data <- object@Data@X[[1]]
  sigma_hat <- lavaan::lav_model_implied(object@Model)$cov[[1]]
  sigma_sqrt <- lavaan::lav_matrix_symmetric_sqrt(sigma_hat)
  s_inv_sqrt <- lavaan::lav_matrix_symmetric_sqrt(object@SampleStats@icov[[1]])
  transformed <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% sigma_sqrt)

  colnames(transformed) <- object@Data@ov.names[[1]]
  errors <- 0 # Not in use for the moment.

  replicate(n_reps, {
    object_ = NULL
    while (is.null(object_)) {
      object_ = tryCatch({
        boot_sample <- transformed[sample(nrow(data), replace = T), ]
        object_ <- lavaan::sem(object, boot_sample)
        stopifnot(lavaan::inspect(object_, "converged"))
        object_
      }, error = function(e) {
        errors <<- errors + 1
        NULL
      })
    }
    progress()
    pvalues(object_)
  })

}
