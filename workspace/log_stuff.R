library("semselector")
library("progressr")
handlers(global = TRUE)

set.seed(22)
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 300
object <- lavaan::sem(model, psych::bfi[1:n, 1:10])
data <- object@Data@X[[1]]
sigma_hat <- lavaan::lav_model_implied(object@Model)$cov[[1]]
sigma_sqrt <- lavaan::lav_matrix_symmetric_sqrt(sigma_hat)
s_inv_sqrt <- lavaan::lav_matrix_symmetric_sqrt(object@SampleStats@icov[[1]])
transformed <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% sigma_sqrt)
colnames(transformed) <- object@Data@ov.names[[1]]
boot_sample <- transformed[sample(nrow(data), replace = T), ]
object <- lavaan::sem(object, boot_sample)


df <- lavaan::fitmeasures(object, "df")
ug <- lavaan::inspect(object, "UG")
eigs <- Re(eigen(ug)$values[1:df])
eigs



n_reps <- 1000
f = function(n_reps) {
  progress <- progressr::progressor(n_reps)
  data <- object@Data@X[[1]]
  sigma_hat <- lavaan::lav_model_implied(object@Model)$cov[[1]]
  sigma_sqrt <- lavaan::lav_matrix_symmetric_sqrt(sigma_hat)
  s_inv_sqrt <- lavaan::lav_matrix_symmetric_sqrt(object@SampleStats@icov[[1]])
  transformed <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% sigma_sqrt)

  colnames(transformed) <- object@Data@ov.names[[1]]
  errors <- 0 # Not in use for the moment.

  replicate(n_reps, {
    result <- NULL
    while (is.null(result)) {
      result <- tryCatch({
        boot_sample <- transformed[sample(nrow(data), replace = T), ]
        object_ <- lavaan::sem(object, boot_sample)
        stopifnot(lavaan::inspect(object_, "converged"))
        progress()
        df <- lavaan::fitmeasures(object_, "df")
        ug <- lavaan::inspect(object_, "UG")
        eigs <- Re(eigen(ug)$values[1:df])
        eigs
      },
      error = function(e) {
        errors <<- errors + 1
        NULL
      }
      )
    }
    result
  })
}
f(n_reps) -> eigens


rowMeans(eigens)
eigens_pred = rowMeans(eigens)[1] - log(1:34)
plot(rowMeans(eigens), eigens_pred)


lambda_log_1 <- mean(eigs) + mean(log(seq(df)))
lambdas_log <- lambda_log_1 - log(seq(df))


rowMeans(eigens)
cov(t(eigens))[1, ] * n
plot(1:34, rowMeans(eigens))


cor(t(eigens))[1, ]



y = rowMeans(eigens)
x = 1:34
summary(lm(y ~ x))$r.squared
summary(lm(y ~ I(x^(-1/3))))$r.squared
summary(lm(y ~ 1 + log(x)))$r.squared

eigens_pred = rowMeans(eigens)[1] - log(1:34)

m = length(eigs)
eig = eigs
k <- ceiling(m / 2)
eig[1:k] <- mean(eig[1:k])
eig[(k + 1):m] <- mean(eig[(k + 1):m])
plot(x, y)
lines(x, predict(lm(y ~ 1 + log(x))))
lines(x, lambdas_log, col = "blue")
points(x, eig, col = "purple", pch = 20)
points(x, eigs, col = "red")
