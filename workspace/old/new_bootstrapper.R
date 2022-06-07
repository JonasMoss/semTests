library("progressr")
handlers(global = TRUE)
handlers(list(
  handler_progress(
    format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta",
    width    = 60,
    complete = "+"
  )
))

semselector(object, 5000)


set.seed(313)
boots <- bootstrapper(object, n_reps = 100)

ggplot2::theme_set(ggplot2::theme_minimal())
ggplot2::ggplot(reshape2::melt(as.data.frame(t(boots))), ggplot2::aes(value)) +
  ggplot2::geom_histogram() +
  ggplot2::xlab("p-value") +
  ggplot2::facet_wrap(~variable, nrow = 3)
