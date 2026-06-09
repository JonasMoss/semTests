# Tests for the internal numerical helpers and their defensive / fallback
# branches: the generalized inverse and orthogonal complement (lavaan_helper.R),
# the provenance print footer (pvalues.R), and the FIML guards (fiml_fmg.R).
# These are reached only on edge cases in normal use, but the edge cases are
# real, so we exercise them directly.

hs  <- " visual =~ x1 + x2 + x3
         textual =~ x4 + x5 + x6
         speed =~ x7 + x8 + x9 "
hs0 <- " visual =~ x1 + a*x2 + x3
         textual =~ x4 + a*x5 + a*x6
         speed =~ x7 + a*x8 + x9 "

# --- generalized_inverse ---------------------------------------------------

test_that("generalized_inverse: full rank, rank-deficient, zero, coercion, bad input", {
  set.seed(1)
  A <- matrix(stats::rnorm(9), 3, 3)
  expect_equal(generalized_inverse(A), solve(A), tolerance = 1e-10)   # full-rank path

  # Moore-Penrose property on a rank-deficient matrix (the partial-rank branch).
  R <- outer(c(1, 2, 3), c(1, 1))                                     # rank 1
  expect_equal(R %*% generalized_inverse(R) %*% R, R, tolerance = 1e-8)

  # all-zero matrix: no positive singular values (the !any(positive) branch).
  expect_true(all(generalized_inverse(matrix(0, 3, 3)) == 0))

  # a non-matrix numeric is coerced via as.matrix().
  gv <- generalized_inverse(c(1, 2, 3))
  expect_equal(dim(gv), c(1L, 3L))

  # input validation.
  expect_error(generalized_inverse(matrix("a", 2, 2)), "numeric")
  expect_error(generalized_inverse(array(0, c(2, 2, 2))), "matrix")
})

# --- get_orthogonal_complement ---------------------------------------------

test_that("get_orthogonal_complement is orthogonal to its input with the right rank", {
  set.seed(2)
  M <- matrix(stats::rnorm(8), 4, 2)                 # 4x2, rank 2
  C <- get_orthogonal_complement(M)
  expect_equal(dim(C), c(4L, 2L))                    # n - rank = 4 - 2
  expect_lt(max(abs(crossprod(C, M))), 1e-10)        # t(C) %*% M == 0
  expect_lt(max(abs(crossprod(C) - diag(2))), 1e-10) # columns orthonormal
})

# --- print provenance footer (pvalues.R) -----------------------------------

test_that("the print footer reports estimator/data/df for single and nested fits", {
  expect_output(print(pvalues(object, "PEBA4")), "estimator:")
  expect_output(print(pvalues(object, "PEBA4")), "df:")
  expect_output(print(pvalues_nested(m0_no_groups, m1_no_groups, tests = "PALL")),
                "nested")
})

test_that("the print footer flags FIML and the nested A.method", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(5); HSm$x1[sample(nrow(HS), 50)] <- NA
  f1 <- lavaan::cfa(hs,  HSm, missing = "fiml", estimator = "MLR")
  f0 <- lavaan::cfa(hs0, HSm, missing = "fiml", estimator = "MLR")
  expect_output(print(pvalues(f1, "PEBA4")), "FIML")              # single-model FIML label
  expect_output(print(pvalues_nested(f0, f1)), "A.method")        # nested FIML footer
})

# --- FIML defensive guards (fiml_fmg.R) ------------------------------------

test_that("fiml_lambdas returns an empty spectrum when df <= 0", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(6); HSm$x1[sample(nrow(HS), 50)] <- NA
  fit <- lavaan::cfa(hs, HSm, missing = "fiml", estimator = "MLR")
  expect_length(fiml_lambdas(fit, 0L)$ug_biased, 0L)
})

test_that("nested FIML requires identical raw data and missingness mask", {
  HS  <- lavaan::HolzingerSwineford1939
  HSa <- HS; set.seed(7); HSa$x1[sample(nrow(HS), 40)] <- NA
  HSb <- HS; set.seed(8); HSb$x2[sample(nrow(HS), 40)] <- NA   # different mask
  f1 <- lavaan::cfa(hs,  HSa, missing = "fiml", estimator = "MLR")
  f0 <- lavaan::cfa(hs0, HSb, missing = "fiml", estimator = "MLR")
  expect_error(pvalues_nested(f0, f1), "same raw data")
})
