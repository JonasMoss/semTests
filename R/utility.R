#' Split vector `x` into `n` chunks of equal size.
#' @param x Input vector.
#' @param n Desired output length.
#' @return List of `n` vectors.
#' @keywords internal
chunk <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

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
