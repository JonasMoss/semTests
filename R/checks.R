# Defensive checks on the fitted object, kept separate from the p-value engine.

#' Warn when a FIML fit uses expected information.
#'
#' Under missing data the expected (Fisher) information is consistent only under
#' MCAR; under MAR -- the regime FIML is adopted for -- the observed information
#' is required for valid inference (Kenward & Molenberghs, 1998). lavaan defaults
#' to observed information for `missing = "ml"`, but it silently accepts
#' `information = "expected"`, which would make the eigenvalue p-values rest on
#' the stronger MCAR assumption. This emits a warning (not an error: expected is
#' still valid under MCAR).
#'
#' lavaan stores `information` as a length-2 vector, so we test its first entry.
#'
#' @param fit A fitted `lavaan` object.
#' @return `fit`, invisibly.
#' @keywords internal
warn_fiml_information <- function(fit) {
  if (!identical(fit@Options$missing, "listwise") &&
      identical(fit@Options$information[1], "expected")) {
    warning(
      "Missing-data (FIML) fit with expected information, which is valid only ",
      "under MCAR. Refit with information = \"observed\" (lavaan's default for ",
      "FIML) for inference valid under MAR.",
      call. = FALSE
    )
  }
  invisible(fit)
}
