#' Split vector `x` into `n` chunks of equal size.
#' @param x Input vector.
#' @param n Desired output length.
#' @return List of `n` vectors.
#' @keywords internal
chunk <- \(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

#' Calculate a saturated model.
#' @keywords internal
#' @param object A `lavaan` object.
#' @return A fitted saturated model.
get_saturated <- \(object) {
  data <- lavaan::lavInspect(object, "data", drop.list.single.group = T)
  data_g <- lapply(seq_along(data), \(x) {
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

#' Returns if not NA; else converts NA to NULL.
#' @keywords internal
nanull <- \(x) {
  if (is.na(x)) {
    NULL
  } else {
    x
  }
}

#' Common default value of 2.
#' @keywords internal
default <- \(x) {
  if (x != "") as.numeric(x) else 2
}


#' Split string into options.
#' @param string Input string
#' @keywords internal

split_input <- \(string) {
  string <- tolower(string)
  splitted <- strsplit(string, "_")[[1]]
  trad <- peba <- eba <- pols <- NULL
  unbiased <- 1
  chisq <- "ml"
  type <- nanull(splitted[1])

  if (length(splitted) == 3) {
    unbiased <- if (splitted[2] == "ug") 2 else 1
    chisq <- splitted[3]
  } else if (length(splitted) == 2) {
    if (splitted[2] == "rls" || splitted[2] == "ml") {
      chisq <- splitted[2]
    } else if (splitted[2] == "ug") {
      unbiased <- 2
    }
  }

  if (startsWith(type, "peba")) {
    peba <- default(substring(type, 5))
  } else if (startsWith(type, "eba")) {
    eba <- default(substring(type, 4))
  } else if (startsWith(type, "pols")) {
    pols <- default(substring(type, 5))
  } else if (type %in% c("std", "sb", "sf", "ss")) {
    trad <- type
  } else {
    stop("Invalid input string in `test`.")
  }

  list(
    trad = trad,
    eba = eba,
    peba = peba,
    pols = pols,
    unbiased = unbiased,
    chisq = chisq
  )
}
