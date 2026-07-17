# Small condition helpers for failures that form part of the public API.
#
# The package still preserves informative errors raised by lavaan and numerical
# internals. Entry-point validation uses stable subclasses so downstream code
# can distinguish bad input, incompatible fits, and fit-quality problems
# without matching English text.

#' Signal a typed semTests error.
#' @keywords internal
semtests_abort <- function(message, subclass = "semTests_error") {
  stop(structure(
    list(message = message, call = NULL),
    class = unique(c(subclass, "semTests_error", "error", "condition"))
  ))
}

#' Signal a typed semTests warning.
#' @keywords internal
semtests_warn <- function(message, subclass = "semTests_warning") {
  warning(structure(
    list(message = message, call = NULL),
    class = unique(c(subclass, "semTests_warning", "warning", "condition"))
  ))
}
