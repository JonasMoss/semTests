### ============================================================================
### Fit model on a subset of the bfi dataset.
### ============================================================================

model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 200
object <- lavaan::sem(model, psych::bfi[1:n, 1:10])
pvalues(object)

### ============================================================================
### Bootstrap setup.
### ============================================================================

# Transform sample for bootstrapping draws
t <- psych::bfi[1:n, 1:10]
f <- lavaan::sem(model, psych::bfi[1:n, 1:10])
sigma_hat <- lavaan::lav_model_implied(object@Model)$cov[[1]]
sigma_sqrt <- lavaan::lav_matrix_symmetric_sqrt(sigma_hat[[1]])
s_inv_sqrt <- lavaan::lav_matrix_symmetric_sqrt(f@SampleStats@icov[[1]])
transformed <- data.frame(as.matrix(t) %*% s_inv_sqrt %*% sigma_sqrt)
colnames(transformed) <- colnames(t)

# Bootstrap setup.
set.seed(1)

# Size of simulation.
n_reps <- 50

# Preparing vectors for data
boots = matrix(NA, nrow = n_reps, ncol = 6)

# Keeping track of the number of errors.
errors <- 0

for (i in seq(n_reps)) {

  if (i %% 50 == 0) print(i)
  result = NULL
  while (is.null(result)) {
    result = tryCatch({
      boot_sample <- transformed[sample(n, replace = T), ]
      model <- lavaan::sem(model, boot_sample)
      stopifnot(lavaan::inspect(model, "converged"))
      list(model = model, boot_sample = boot_sample)
    }, error = function(e) {
      errors <<- errors + 1
      NULL
    })
  }
  model = result$model
  boot_sample = result$boot_sample

  boots[i, ] = pvalues(model)
}


### ============================================================================
### Plotting and saving
### ============================================================================

#saveRDS(my.df, "illustration1_R5000.rds")

ggplot2::theme_set(ggplot2::theme_minimal())
ggplot2::ggplot(reshape2::melt(pvalues), ggplot2::aes(value)) +
  ggplot2::geom_histogram() +
  ggplot2::xlab("p-value") +
  ggplot2::facet_wrap(~variable, nrow = 3)

# Distances
kolmogorov_smirnov <- as.vector(sapply(pvalues, function(x) ks.test(x, y = "punif")$statistic))

max(abs(sort(x) - (0:99)/100))

anderson_darling <- as.vector(sapply(pvalues, function(x) goftest::ad.test(x, null = "punif")$statistic))

cat("Kolmogorov-Smirnov test selects: ", colnames(pvalues)[which.min(anderson_darling)], "\n")
cat("Anderson-Darling test selects: ", colnames(pvalues)[which.min(kolmogorov_smirnov)])
