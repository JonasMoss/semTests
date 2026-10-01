# The nested restriction matrix `A` (Satorra, 2000). Its rows span the
# directions of H1's effective parameter space that H0 removes, so it is the
# orthogonal complement of the tangent map J from H0's effective parameters
# into H1's. Three constructions, chosen by `A.method`:
#
# "delta"  (default) J = Delta_1^+ Delta_0 with both Jacobians at the H0
#          solution: Delta_1 is evaluated at the H1 parameter vector that
#          reproduces H0's fitted moments. When H0 is a parameter map
#          theta_1 = g(theta_0), the chain rule gives
#          Delta_0(theta_0) = Delta_1(g(theta_0)) J, so this is J itself and A
#          is exact. It needs no parameter correspondence beyond locating the
#          H0 solution in H1.
# "exact"  J from the two models' parameter bases aligned by parameter keys.
#          Requires the same underlying full free-parameter set.
# "legacy" J = Delta_1^+ Delta_0 with Delta_1 at the H1 estimate and Delta_0
#          at the H0 estimate, as in semTests <= 1.0.0 and lavaan's default
#          lavTestLRT(A.method = "delta"). The Jacobians are taken at points
#          O(n^-1/2) apart, so A, and the reference spectrum, carry an
#          O(n^-1/2) error that the other two avoid.
#
# lavaan exposes Delta directly via lavInspect(fit, "delta"), so we read it
# from there (for "delta", on a copy of H1 moved to the H0 solution).

#' Nested restriction matrix in H1's effective parameter space.
#' @param m1,m0 The unrestricted and restricted `lavaan` fits.
#' @param method `"delta"`, `"exact"`, or `"legacy"`; see above.
#' @keywords internal
get_a_matrix <- function(m1, m0, method = c("delta", "exact", "legacy")) {
  method <- match.arg(method)
  if (method == "exact") {
    basis1 <- model_parameter_basis(m1)
    basis0 <- align_parameter_basis(m0, m1, model_parameter_basis(m0))
    return(t(get_orthogonal_complement(generalized_inverse(basis1) %*% basis0)))
  }
  at <- m1
  if (method == "delta") {
    at <- h1_at_h0_solution(m1, m0)
    if (is.null(at)) {
      semtests_warn(
        paste0(
          "Could not place the restricted solution in the unrestricted ",
          "model's parameter space, so the restriction matrix falls back to ",
          "`A.method = \"legacy\"`, whose reference spectrum carries an ",
          "O(n^-1/2) error. See `?pvalues_nested`."
        ),
        "semTests_warning_legacy_restriction"
      )
      at <- m1
    }
  }
  map <- generalized_inverse(model_effective_delta(at)) %*%
    model_effective_delta(m0)
  t(get_orthogonal_complement(map))
}

#' The unrestricted fit moved to the restricted solution.
#'
#' Locates the H1 parameter vector that reproduces H0's fitted moments by
#' matching parameter keys: each free H1 parameter takes H0's value of the
#' same parameter, and a parameter H0 omits entirely is zero. The candidate
#' is accepted only if H1's implied moments there equal H0's fitted moments.
#' Returns a copy of `m1` at that point, or `NULL` when no such point is found.
#' @keywords internal
h1_at_h0_solution <- function(m1, m0, tolerance = 1e-6) {
  keys <- function(pt) {
    fields <- intersect(c("lhs", "op", "rhs", "group", "level", "block"), names(pt))
    do.call(paste, c(unname(pt[fields]), sep = "\r"))
  }
  pt1 <- lavaan::parTable(m1)
  pt0 <- lavaan::parTable(m0)
  free1 <- pt1[pt1$free > 0L, , drop = FALSE]
  x <- lavaan::lav_model_get_parameters(m1@Model)
  if (!nrow(free1) || length(x) != m1@Model@nx.free || max(free1$free) > length(x)) {
    return(NULL)
  }
  hit <- match(keys(free1), keys(pt0))
  value <- pt0$est[hit]
  value[is.na(hit)] <- 0
  x[free1$free] <- value
  moved <- m1
  moved@Model <- lavaan::lav_model_set_parameters(m1@Model, x = x)
  implied1 <- lavaan::lav_model_implied(moved@Model)
  implied0 <- lavaan::lav_model_implied(m0@Model)
  for (part in intersect(names(implied0), names(implied1))) {
    a <- unlist(implied0[[part]])
    b <- unlist(implied1[[part]])
    if (!is.numeric(a)) next
    if (length(a) != length(b) || anyNA(b) ||
      max(abs(a - b)) > tolerance * max(1, abs(a))) {
      return(NULL)
    }
  }
  moved
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
