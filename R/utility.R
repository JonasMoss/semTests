
#' Is this the classical normal-theory, complete-data, ML case?
#'
#' The Du-Bentler unbiased gamma and the RLS (`browne.residual.nt.model`)
#' statistic are only defined here; off this case the eigenvalue machinery uses
#' the biased gamma and the estimator's own (uncorrected) statistic.
#' @keywords internal
is_classic_nt <- function(fit) {
  fit@Options$estimator == "ML" &&
    !isTRUE(fit@Model@categorical) &&
    identical(fit@Options$missing, "listwise")
}

#' @keywords internal
make_chisqs <- function(chisq, m0, m1) {
  # lavaan >= 0.7-1 returns a named list of tests from `lavTest()` (e.g. both
  # "standard" and the requested test), whereas earlier versions return a single
  # flat test object with `$stat`. Grab the statistic from whichever layout we get.
  get_stat <- function(object, test) {
    res <- lavaan::lavTest(object, test = test)
    if (!is.null(res[[test]])) res[[test]][["stat"]] else res[["stat"]]
  }
  # "ml" is the standard (uncorrected) test statistic of the fit -- the LRT for
  # ML, and the corresponding uncorrected quadratic form for ULS/GLS/DWLS/FIML.
  # It is the right input to the eigenvalue correction for any estimator; the
  # native scaled/shifted statistic is already corrected and must NOT be used.
  standard <- function(object) get_stat(object, "standard")
  rls <- function(object) get_stat(object, "browne.residual.nt.model")
  wrap <- function(f, object) if (missing(object)) 0 else f(object)

  classic_nt <- is_classic_nt(m0)
  # The default statistic is fit-appropriate: RLS for the classical case (its
  # historical default), the standard statistic otherwise.
  chisq[chisq == "auto"] <- if (classic_nt) "rls" else "ml"
  chisq <- unique(chisq)

  if ("rls" %in% chisq && !classic_nt) {
    stop("The RLS statistic ('browne.residual.nt.model') is only available for ",
         "continuous, complete-data ML; it silently degrades to ADF otherwise. ",
         "Use the standard statistic (the `_ML` suffix) or omit the suffix.",
         call. = FALSE)
  }

  chisqs <- c()
  if ("ml" %in% chisq) chisqs["ml"] <- standard(m0) - wrap(standard, m1)
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
  chisq <- "auto"   # resolved per fit in make_chisqs (rls if classical NT, else ml)
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
