# The nested restriction matrix `A` (Satorra, 2000). Its rows span the
# directions of H1's effective parameter space that H0 removes, so it is the
# orthogonal complement of the tangent map J from H0's effective parameters
# into H1's. Two constructions, chosen by `A.method` and named as in lavaan's
# lavTestLRT():
#
# "delta"  (default) J = Delta_1^+ Delta_0 with each model's Jacobian at its
#          own estimate, lavaan's default. Built from moment-space Jacobians,
#          it applies to any nesting, including H0s that drop or fix paths,
#          and does not depend on how the factors are identified.
# "exact"  J from the two models' parameter bases aligned by parameter keys,
#          lavaan's "exact". Requires the same underlying full free-parameter
#          set in both models. The map is exact in one parameterization while
#          the other ingredients sit at the H1 estimate, so the reference
#          spectrum moves at O(n^-1/2) with the identification (std.lv versus
#          marker loadings).
#
# lavaan exposes Delta directly via lavInspect(fit, "delta"), so we read it
# from there instead of re-deriving it.

#' Nested restriction matrix in H1's effective parameter space.
#' @param m1,m0 The unrestricted and restricted `lavaan` fits.
#' @param method `"delta"` or `"exact"`; see above.
#' @keywords internal
get_a_matrix <- function(m1, m0, method = c("delta", "exact")) {
  method <- match.arg(method)
  if (method == "exact") {
    basis1 <- model_parameter_basis(m1)
    basis0 <- align_parameter_basis(m0, m1, model_parameter_basis(m0))
    return(t(get_orthogonal_complement(generalized_inverse(basis1) %*% basis0)))
  }
  map <- generalized_inverse(model_effective_delta(m1)) %*%
    model_effective_delta(m0)
  t(get_orthogonal_complement(map))
}

#' Align an H0 parameter basis to H1's full-parameter ordering.
#' @keywords internal
align_parameter_basis <- function(m0, m1, basis0) {
  key0 <- full_parameter_keys(m0)
  key1 <- full_parameter_keys(m1)
  if (length(key0) != length(key1) ||
    anyDuplicated(key0) || anyDuplicated(key1) ||
    !setequal(key0, key1)) {
    stop("The exact nested restriction requires the same underlying full ",
      "free-parameter set in both models. Use `A.method = \"delta\"` for ",
      "models represented by different parameter sets.",
      call. = FALSE
    )
  }
  basis0[match(key1, key0), , drop = FALSE]
}

#' Keys of lavaan's full free parameters, in free-parameter order.
#' @keywords internal
full_parameter_keys <- function(fit) {
  pt <- lavaan::parTable(fit)
  pt <- pt[pt$free > 0L, , drop = FALSE]
  pt <- pt[order(pt$free), , drop = FALSE]
  if (nrow(pt) != fit@Model@nx.free ||
    !identical(as.integer(pt$free), seq_len(fit@Model@nx.free))) {
    stop("Could not identify lavaan's full free-parameter ordering.",
      call. = FALSE
    )
  }
  fields <- intersect(
    c("lhs", "op", "rhs", "group", "level", "block"),
    names(pt)
  )
  do.call(paste, c(unname(pt[fields]), sep = "\r"))
}
