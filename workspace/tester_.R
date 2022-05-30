library("semselector")
set.seed(313)
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 200
object <- lavaan::sem(model, psych::bfi[1:n, 1:10])
selector <- semselector(object, n_reps = 50)
options(warn=-1)





data <- bollen_stine_transform(object)
boot_object <- bootstrap(object, data[[1]])

pvalues_one(boot_object)
