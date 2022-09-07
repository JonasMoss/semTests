test_that("distance function works", {
  x <- seq(0.001, 0.999, length.out = 100)
  dist <- c(
    "kolmogorov-smirnov",
    "anderson-darling",
    "kullback-leibler",
    "cramer-von mises",
    "0.05-distance"
  )
  expect_equal(length(sapply(dist, function(d) distance(x, d))), length(dist))
})
