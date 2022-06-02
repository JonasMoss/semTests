test_that("bootstrapper runs, two models.", {
  boots = bootstrapper(m0, m1, functional = identity, n_reps = 1)
  expect_equal(dim(boots), c(2, 1))
})

test_that("bootstrapper runs, one models", {
  boots = bootstrapper(m0, m1 = NULL, functional = identity, n_reps = 2)
  expect_equal(length(boots), 2)
})


test_that("bollens_stine_transform works", {
  s <- s_and_s_inv(object)

  lhs <- lapply(seq(object@SampleStats@ngroups), function(i) {
    data <- object@Data@X[[i]]
    s_sqrt <- s[[i]]$s_sqrt
    s_inv_sqrt <- s[[i]]$s_inv_sqrt
    frame <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% s_sqrt)
    colnames(frame) <- object@Data@ov.names[[i]]
    frame
  })

  rhs <- bollen_stine_transform(object)
  expect_equal(lhs, rhs)
})

test_that("bootstrap runs", {
  set.seed(313)
  boot_1 = bootstrap(m0, m0, data)
  set.seed(313)
  boot_2 = bootstrap(m0, data = data)
  expect_equal(boot_2, boot_2)
})

test_that("s_and_s_inv works", {
  lhs <- lapply(seq(object@SampleStats@ngroups), function(i) {
    s_hat <- lavaan::lav_model_implied(object@Model)$cov[[i]]
    s_inv_hat <- object@SampleStats@icov[[i]]
    list(
      s_sqrt = lavaan::lav_matrix_symmetric_sqrt(s_hat),
      s_inv_sqrt = lavaan::lav_matrix_symmetric_sqrt(s_inv_hat)
    )
  })

  rhs <- s_and_s_inv(object)
  expect_equal(lhs, rhs)
})
