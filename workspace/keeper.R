
model <- "A =~ A1+A2+A3+A4+A5;
          C =~ C1+C2+C3+C4+C5"
n <- 50

#lavaan::lavaan(model = model, data = psych::bfi[1:n, 1:10])


m0 <- lavaan::cfa(model, psych::bfi[1:n, 1:10])

data <- bollen_stine_transform(m0)

ns <- m0@Data@nobs
ids <- lapply(ns, function(n) sample(x = n, size = n, replace = TRUE))

boot_sample <- lavaan::lav_data_update(
  lavdata = m0@Data,
  newX = lapply(seq_along(ns), function(i) data[[i]][ids[[i]], ]),
  lavoptions = lavaan::lavInspect(m0, "options")
)

boot_m0 <- lavaan::lavaan(
  slotOptions = m0@Options,
  slotParTable = m0@ParTable,
  slotData = boot_sample,
  control = list(iter.max = 10)
)


boot_m0
