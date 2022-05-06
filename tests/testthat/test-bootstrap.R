test_that("bollens_stine_transform works", {
  models <- list(m0)
  s <- s_and_s_inv(models[[1]])
  lhs <- lapply(models, function(model) {
    lapply(seq(model@SampleStats@ngroups), function(i) {
      data <- model@Data@X[[i]]
      s_sqrt <- s[[i]]$s_sqrt
      s_inv_sqrt <- s[[i]]$s_inv_sqrt
      frame <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% s_sqrt)
      colnames(frame) <- model@Data@ov.names[[i]]
      frame
    })
  })
  rhs <- bollen_stine_transform(m0)
  expect_equal(lhs, rhs)

  models <- list(m0, m1)
  s <- s_and_s_inv(models[[1]])
  lhs <- lapply(models, function(model) {
    lapply(seq(model@SampleStats@ngroups), function(i) {
      data <- model@Data@X[[i]]
      s_sqrt <- s[[i]]$s_sqrt
      s_inv_sqrt <- s[[i]]$s_inv_sqrt
      frame <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% s_sqrt)
      colnames(frame) <- model@Data@ov.names[[i]]
      frame
    })
  })
  rhs <- bollen_stine_transform(m0, m1)
  expect_equal(lhs, rhs)
})

test_that("bootstrap works", {
  set.seed(1)
  data <- object@Data@X
  ns <- object@Data@nobs
  ids <- lapply(ns, function(n) sample(x = n, size = n, replace = TRUE))

  boot_sample <- lavaan::lav_data_update(
    lavdata = object@Data,
    newX = lapply(seq_along(ns), function(i) data[[i]][ids[[i]], ]),
    lavoptions = lavaan::lavInspect(object, "options")
  )

  boot_object <- lavaan::lavaan(
    slotOptions = object@Options,
    slotParTable = object@ParTable,
    slotData = boot_sample
  )
  set.seed(1)
  expect_equal(boot_object@ParTable, bootstrap(object, data)@ParTable)
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
