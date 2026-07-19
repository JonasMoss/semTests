# The nested restriction matrix `A` is the orthogonal complement of the map
# between the two models' Jacobians (Delta). lavaan exposes Delta directly via
# lavInspect(fit, "delta") -- bit-for-bit the matrix the formerly-vendored
# computeDelta()/derivative.*.LISREL machinery recomputed -- so we read it from
# there instead of re-deriving it. (The single-model FIML path likewise reads
# Delta from lavInspect; see fiml_fmg.R.)

get_a_matrix <- function(m1, m0) {
  map <- generalized_inverse(model_effective_delta(m1)) %*%
    model_effective_delta(m0)
  t(get_orthogonal_complement(map))
}
