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
  expect_equal(generalized_inverse(A), solve(A), tolerance = 1e-8)   # full-rank path

  # Moore-Penrose property on a rank-deficient matrix (the partial-rank branch).
  R <- outer(c(1, 2, 3), c(1, 1))                                     # rank 1
  expect_equal(R %*% generalized_inverse(R) %*% R, R, tolerance = 1e-8)

  # all-zero matrix: no positive singular values (the !any(positive) branch).
  expect_true(all(generalized_inverse(matrix(0, 3, 3)) == 0))

  # a non-matrix numeric is coerced via as.matrix().
  gv <- generalized_inverse(c(1, 2, 3))
  expect_equal(dim(gv), c(1L, 3L))

  complex_matrix <- matrix(c(1 + 1i, 2 - 1i, 3, 4 + 2i), 2L)
  expect_equal(
    generalized_inverse(complex_matrix),
    solve(complex_matrix),
    tolerance = 1e-8
  )

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

  plain_ml <- lavaan::cfa(
    hs, lavaan::HolzingerSwineford1939,
    estimator = "ML"
  )
  expect_output(print(pvalues(plain_ml, "PEBA4_ML")), "estimator: ML \\|")
})

test_that("the print footer flags FIML and the nested A.method", {
  HS  <- lavaan::HolzingerSwineford1939
  HSm <- HS; set.seed(5); HSm$x1[sample(nrow(HS), 50)] <- NA
  f1 <- lavaan::cfa(hs,  HSm, missing = "fiml", estimator = "MLR")
  f0 <- lavaan::cfa(hs0, HSm, missing = "fiml", estimator = "MLR")
  expect_output(print(pvalues(f1, "PEBA4")), "FIML")              # single-model FIML label
  expect_output(print(pvalues_nested(f0, f1)), "A.method")        # nested FIML footer
})

test_that("estimator provenance falls back cleanly when estimator.orig is absent", {
  fit <- object
  fit@Options$estimator.orig <- NULL
  expect_equal(requested_estimator(fit), fit@Options$estimator)
  expect_equal(
    fit_provenance(fit, nested = FALSE)$estimator_requested,
    fit@Options$estimator
  )
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

test_that("FIML helper guards explain malformed fit and matrix shapes", {
  HS <- lavaan::HolzingerSwineford1939
  HSm <- HS
  set.seed(9)
  HSm$x1[sample(nrow(HS), 40)] <- NA
  fit <- lavaan::cfa(hs, HSm, missing = "fiml", estimator = "MLR")

  categorical_data <- HS
  categorical_data[paste0("x", 1:9)] <- lapply(
    categorical_data[paste0("x", 1:9)],
    function(x) ordered(cut(x, 3))
  )
  categorical <- lavaan::cfa(
    hs, categorical_data,
    ordered = paste0("x", 1:9)
  )
  expect_error(
    fiml_check_supported(categorical),
    "only implemented for continuous lavaan fits"
  )
  fixed_x <- lavaan::cfa(
    "visual =~ x1 + x2 + x3
     visual ~ ageyr",
    HSm, missing = "fiml", estimator = "MLR", fixed.x = TRUE
  )
  expect_error(
    fiml_check_supported(fixed_x),
    "fixed or conditional observed exogenous"
  )
  expect_error(
    fiml_group_matrices(list(diag(2), diag(2)), fit, "test"),
    "Expected 1 test matrix"
  )
  invalid_weights <- fit
  invalid_weights@SampleStats@ntotal <- 0L
  expect_error(
    fiml_group_weights(invalid_weights),
    "valid FIML group weights"
  )
  local({
    testthat::local_mocked_bindings(
      fiml_group_matrices = function(...) {
        list(matrix(
          0, nrow = 2L, ncol = 8L,
          dimnames = list(NULL, paste0("x", seq_len(8L)))
        ))
      },
      .package = "semTests"
    )
    expect_error(
      fiml_data_blocks(fit),
      "do not contain all observed variables"
    )
  })
  expect_error(
    fiml_validate_A(matrix(0, 2L, 3L), df = 1L, npar = 3L),
    "rank.*does not match"
  )
  expect_error(
    fiml_validate_A(matrix(0, 1L, 2L), df = 1L, npar = 3L),
    "does not match H1"
  )

  malformed <- fit
  malformed@Model@nx.free <- malformed@Model@nx.free + 1L
  expect_error(
    fiml_parameter_keys(malformed),
    "Could not identify lavaan's full free-parameter ordering"
  )
})

test_that("FIML data and dimension invariants fail close to their source", {
  HSa <- lavaan::HolzingerSwineford1939
  HSb <- HSa
  set.seed(10)
  HSa$x1[sample(nrow(HSa), 30L)] <- NA
  set.seed(11)
  HSb$x2[sample(nrow(HSb), 30L)] <- NA
  f1 <- lavaan::cfa(hs, HSa, missing = "fiml", estimator = "MLR")
  f0 <- lavaan::cfa(hs0, HSb, missing = "fiml", estimator = "MLR")
  expect_error(fiml_check_same_data(f0, f1), "same raw data")

  testthat::local_mocked_bindings(
    fiml_h1_information_observed = function(...) diag(2L),
    fiml_gamma_matrix = function(...) diag(3L),
    fiml_delta_effective = function(...) matrix(0, 2L, 1L),
    .package = "semTests"
  )
  expect_error(
    fiml_lambdas(f1, 1L, fiml.convention = "observed"),
    "dimensions do not agree"
  )
})

test_that("nested FIML checks restriction and information dimensions", {
  HS <- lavaan::HolzingerSwineford1939
  set.seed(12)
  HS$x1[sample(nrow(HS), 30L)] <- NA
  f1 <- lavaan::cfa(hs, HS, missing = "fiml", estimator = "MLR")
  f0 <- lavaan::cfa(hs0, HS, missing = "fiml", estimator = "MLR")

  testthat::local_mocked_bindings(
    fiml_A_delta = function(...) matrix(0, 1L, 1L),
    .package = "semTests"
  )
  expect_error(
    fiml_lambdas_nested(f0, f1, 1L, A.method = "delta"),
    "restriction matrix and H1 information dimensions"
  )
})

test_that("FIML keeps compatibility with lavaan simple-constraint bases", {
  HS <- lavaan::HolzingerSwineford1939
  set.seed(13)
  HS$x1[sample(nrow(HS), 30L)] <- NA
  fit <- lavaan::cfa(hs, HS, missing = "fiml", estimator = "MLR")
  npar <- fit@Model@nx.free

  legacy <- fit
  legacy@Model@eq.constraints <- FALSE
  legacy@Model@ceq.simple.only <- TRUE
  legacy@Model@ceq.simple.K <- diag(npar)[, -npar, drop = FALSE]
  expect_equal(dim(fiml_K_matrix(legacy)), c(npar, npar - 1L))

  empty <- fit
  empty@Model@eq.constraints <- FALSE
  empty@Model@ceq.simple.only <- TRUE
  empty@Model@ceq.simple.K <- matrix(numeric(), 0L, npar)
  expect_equal(fiml_K_matrix(empty), diag(npar))
})
