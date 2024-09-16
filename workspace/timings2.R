# n <- 1000
# A <- matrix(runif(n^2)*2-1, ncol=n)
# Sigma <- t(A) %*% A
#
# k <- 10
#
# microbenchmark::microbenchmark(
#   RSpectra::eigs(Sigma, k = k, which = "LR")$values,
#   RSpectra::eigs(Sigma, k = k, which = "LR", opts = list(retvec = FALSE))$values,
#   eigen(Sigma, only.values = TRUE)$values[seq(k)],
#   eigen(Sigma)$values[seq(k)])


Sigma <- ugamma_nested(m0, m1)[[1]]
Sigma1 <- Sigma
Sigma1[abs(Sigma) < 10^(-9)] = 0
k <- 6
Sigma2 <- Sigma1

#' Create sparse matrix
#' @param mat Matrix input.
#' @param lim Elements with absolute value less than `lim` get set to `0`.
#' @return Object of `dgCMatrix`.
sparsify <- \(mat, lim = 1e-9) {
  mat[abs(mat) < lim] = 0
  Matrix::Matrix(mat, sparse = TRUE)
}

Sigma1[abs(Sigma) < 10^(-9)] = 0
Sigma1 <- Matrix::Matrix(Sigma1, sparse = TRUE)
Sigma2 <- Matrix::Matrix(Sigma)

k <- dim(Sigma1)[1]
k <- 6

microbenchmark::microbenchmark(
  RSpectra::eigs(Sigma1, k = k, which = "LR")$values,
  RSpectra::eigs(Sigma2, k = k, which = "LR")$values,
  RSpectra::eigs(Sigma, k = k, which = "LR")$values,
  eigen(Sigma, only.values = TRUE)$values,
  eigen(Sigma)$values)

microbenchmark::microbenchmark(ugamma_nested(m0, m1)[[1]],
                               eigen(Sigma, only.values = TRUE)$values)


mod = paste0("f=~", paste0("x", 1:30,collapse="+"))
data <- lavaan::simulateData(mod)
data$g <- rep(1:2, each = 250)

m0 <- lavaan::cfa(model = mod, data = data, estimator = "MLM", group = "g", group.equal = "loadings")
m1 <- lavaan::cfa(model = mod, data = data, estimator = "MLM", group = "g")

k <- lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")

ug <- ugamma_nested(m0, m1)[[1]]
ug_sparse <- sparsify(ug)

microbenchmark::microbenchmark(ugamma_nested(m0, m1)[[1]],
                               eigen(ug, only.values = TRUE)$values,
                               RSpectra::eigs(ug_sparse, k = k, which = "LR")$values,
                               times = 10)

for(i in 1:1000) {
  RSpectra::eigs(Sigma, k = k, which = "LR")$values
}


microbenchmark::microbenchmark(
  RSpectra::eigs(Sigma, k = k, which = "LR")$values,
  RSpectra::eigs(Sigma1, k = k, which = "LR")$values)



microbenchmark::microbenchmark(
  RSpectra::eigs(Sigma, k = k, which = "LR")$values,
  RSpectra::eigs(Sigma, k = k, which = "LR", opts = list(retvec = FALSE))$values,
  eigen(Sigma, only.values = TRUE)$values[seq(k)],
  eigen(Sigma)$values[seq(k)])


microbenchmark::microbenchmark(
  RSpectra::eigs(Sigma, k = k, which = "LR")$values,
  RSpectra::eigs(Sigma, k = k, which = "LR", opts = list(retvec = FALSE))$values,
  eigen(Sigma, only.values = TRUE)$values[seq(k)],
  eigen(Sigma)$values[seq(k)])
