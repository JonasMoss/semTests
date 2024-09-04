first <- readRDS("workspace/first.Rds")
second <- readRDS("workspace/second.Rds")


res <- lav_ugamma_nested_2000(first, second, get_gamma(first, collapse = TRUE))
#pvalues_nested(first, second)

object <- first
unbiased <- FALSE
collapse <- FALSE

stopifnot(is.logical(unbiased))
lavoptions = lavaan::lavInspect(object, "options")
lavdata = object@Data
lavoptions$gamma.unbiased <- unbiased
gamma_list <- lavaan:::lav_samplestats_from_data(lavdata, lavoptions)@NACOV
if (collapse) {
  #as.matrix(Matrix::bdiag(gamma_list))
  lavaan:::lav_matrix_bdiag(gamma_list)
} else {
  gamma_list
}
