# The nested restriction matrix `A` is the orthogonal complement of the map
# between the two models' Jacobians (Delta). lavaan exposes Delta directly via
# lavInspect(fit, "delta") -- bit-for-bit the matrix the formerly-vendored
# computeDelta()/derivative.*.LISREL machinery recomputed -- so we read it from
# there instead of re-deriving it. (The single-model FIML path likewise reads
# Delta from lavInspect; see fiml_fmg.R.)

get_a_matrix <- function(m1, m0) {

  delta <- function(m) {
    d <- lavaan::lavInspect(m, "delta")
    d <- if (is.list(d)) do.call(rbind, lapply(d, as.matrix)) else as.matrix(d)
    if (m@Model@eq.constraints) {
      return(d %*% m@Model@eq.constraints.K)
    }
    if (methods::.hasSlot(m@Model, "ceq.simple.only") && m@Model@ceq.simple.only) {
      K <- as.matrix(m@Model@ceq.simple.K)
      return(d %*% qr.Q(qr(K)))
    }
    if (methods::.hasSlot(m@Model, "ceq.JAC") && nrow(m@Model@ceq.JAC) > 0L) {
      K <- get_orthogonal_complement(t(as.matrix(m@Model@ceq.JAC)))
      return(d %*% K)
    }
    d
  }

  t(get_orthogonal_complement(generalized_inverse(delta(m1)) %*% delta(m0)))
}
