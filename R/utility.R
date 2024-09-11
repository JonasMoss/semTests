
#' @keywords internal
make_chisqs <- \(chisq, m0, m1) {
  ml <- \(object) lavaan::lavTest(object, test = "standard")$stat
  rls <- \(object) lavaan::lavTest(object, test = "browne.residual.nt.model")$stat
  wrap <- \(f, object) if (missing(object)) 0 else f(object)
  chisqs <- c()
  if ("ml" %in% chisq) chisqs["ml"] <- ml(m0) - wrap(ml, m1)
  if ("rls" %in% chisq) chisqs["rls"] <- rls(m0) - wrap(rls, m1)
  chisqs
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

#' Create sparse matrix
#' @param mat Matrix input.
#' @param lim Elements with absolute value less than `lim` get set to `0`.
#' @return Object of `dgCMatrix`.
sparsify <- \(mat, lim = 1e-9) {
  mat[abs(mat) < lim] = 0
  Matrix::Matrix(mat, sparse = TRUE)
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
