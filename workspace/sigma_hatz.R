library("semselector")
set.seed(313)
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5;
          E =~ E1+E2+E3+E4+E5;
          N =~ N1+N1+N3+N4+N5"

model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5;"
n <- 200
object <- lavaan::sem(model, psych::bfi[1:n, 1:20])
selector <- semselector(object, n_reps = 5000)
selector
plot(selector, binwidth = 0.05)


data <- object@Data@X[[1]]
sigma_hat <- lavaan::lav_model_implied(object@Model)$cov
sigma_sqrt <- lapply(sigma_hat, lavaan::lav_matrix_symmetric_sqrt)
s_inv_sqrt <- lapply(object@SampleStats@icov, lavaan::lav_matrix_symmetric_sqrt)
transformed <- data.frame(as.matrix(data) %*% s_inv_sqrt %*% sigma_sqrt)
