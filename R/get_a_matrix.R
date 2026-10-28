# The nested restriction matrix `A` (Satorra, 2000). Its rows span the
# directions of H1's effective parameter space that H0 removes, so it is the
# orthogonal complement of the tangent map J from H0's effective parameters
# into H1's. Three constructions, chosen by `A.method`:
#
# "weighted" (default) J = (Delta_1' W Delta_1)^+ Delta_1' W Delta_0, the
#          least-squares map from Delta_0 onto Delta_1 in the fit's own
#          moment-space weight W (the group-weighted WLS.V for complete data,
#          the observed saturated information for FIML), with each model's
#          Jacobian at its own estimate. Delta_0 lies in the span of Delta_1
#          only at a common point, so the projection metric matters: W
#          transforms inversely to Delta under a change of units, which makes
#          the test invariant to both the units of the observed variables and
#          the identification of the factors.
# "delta"  J = Delta_1^+ Delta_0, the same map in the Euclidean metric. This is
#          lavaan's default in lavTestLRT(method = "satorra.2000") and
#          reproduces it. It does not depend on the identification, but it
#          does depend on the units of the observed variables.
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
#' @param method `"weighted"`, `"delta"`, or `"exact"`; see above.
#' @keywords internal
get_a_matrix <- function(m1, m0, method = c("weighted", "delta", "exact")) {
  method <- match.arg(method)
  if (method == "exact") {
    basis1 <- model_parameter_basis(m1)
    basis0 <- align_parameter_basis(m0, m1, model_parameter_basis(m0))
    return(t(get_orthogonal_complement(generalized_inverse(basis1) %*% basis0)))
  }
  delta1 <- model_effective_delta(m1)
  delta0 <- model_effective_delta(m0)
  if (method == "delta") {
    return(t(get_orthogonal_complement(generalized_inverse(delta1) %*% delta0)))
  }
  # Work in standardized parameter coordinates, where Delta' W Delta has a
  # unit diagonal: the rank decisions below then do not depend on the units of
  # the parameters. A covector a on the standardized H1 parameters is a * s1 on
  # the original ones, and the test does not depend on the scale of A's rows.
  weight <- restriction_weight(m1)
  s1 <- positive_scale(colSums(delta1 * (weight %*% delta1)))
  s0 <- positive_scale(colSums(delta0 * (weight %*% delta0)))
  delta1 <- sweep(delta1, 2L, s1, "/")
  weighted <- weight %*% delta1
  map <- generalized_inverse(crossprod(delta1, weighted)) %*%
    crossprod(weighted, sweep(delta0, 2L, s0, "/"))
  sweep(t(get_orthogonal_complement(map)), 2L, s1, "*")
}

#' Moment-space weight of the weighted restriction map.
#'
#' The group-weighted `WLS.V` for complete data, as in [nested_factor_2000()],
#' and the observed saturated information for FIML, as in the FIML nested
#' spectrum.
#' @keywords internal
restriction_weight <- function(fit) {
  if (is_fiml(fit)) return(fiml_h1_information_observed(fit))
  weights <- lavaan::lavTech(fit, "WLS.V")
  shares <- unlist(fit@SampleStats@nobs) / fit@SampleStats@ntotal
  lavaan::lav_matrix_bdiag(Map(`*`, shares, weights))
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
      "free-parameter set in both models. Use `A.method = \"weighted\"` for ",
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
