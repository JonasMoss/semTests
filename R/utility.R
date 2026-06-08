
#' @keywords internal
make_chisqs <- function(chisq, m0, m1) {
  # lavaan >= 0.7-1 returns a named list of tests from `lavTest()` (e.g. both
  # "standard" and the requested test), whereas earlier versions return a single
  # flat test object with `$stat`. Grab the statistic from whichever layout we get.
  get_stat <- function(object, test) {
    res <- lavaan::lavTest(object, test = test)
    if (!is.null(res[[test]])) res[[test]][["stat"]] else res[["stat"]]
  }
  ml <- function(object) get_stat(object, "standard")
  rls <- function(object) get_stat(object, "browne.residual.nt.model")
  wrap <- function(f, object) if (missing(object)) 0 else f(object)
  chisqs <- c()
  if ("ml" %in% chisq) chisqs["ml"] <- ml(m0) - wrap(ml, m1)
  if ("rls" %in% chisq) chisqs["rls"] <- rls(m0) - wrap(rls, m1)
  chisqs
}

#' Returns if not NA; else converts NA to NULL.
#' @keywords internal
nanull <- function(x) {
  if (is.na(x)) {
    NULL
  } else {
    x
  }
}

#' Common default value of 2.
#' @keywords internal
default <- function(x) {
  if (x != "") as.numeric(x) else 2
}

#' Split string into options.
#' @param string Input string
#' @keywords internal

split_input <- function(string) {

  string <- tolower(string)
  splitted <- strsplit(string, "_")[[1]]
  trad <- peba <- eba <- pols <- NULL
  unbiased <- 1
  chisq <- "rls"
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
  } else if (type %in% c("std", "sb", "sf", "ss", "all", "pall")) {
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
