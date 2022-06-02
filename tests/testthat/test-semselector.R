test_that("semselector works for one model", {
  set.seed(313)
  sel = semselector(m1, n_reps = 2, distances = "anderson-darling")
  expect_equal(dim(sel), c(1, 3))
})
