#' Defunct estimators
#'
#' @description
#' `r lifecycle::badge("defunct")`
#'
#' These estimators were deprecated to encourage use of `lmtp_tmle()` and `lmtp_sdr()` functions.
#'
#' @keywords internal
#' @name defunct
NULL

#' @export
#' @rdname defunct
lmtp_sub <- function(...) {
  lifecycle::deprecate_stop(
    when = "1.5.0",
    what = "lmtp_sub()",
    details = c("G-computation requires the use of correctly pre-specified parametric models for valid statistical inference. Use `lmtp_tmle()` or `lmtp_sdr()`.")
  )
}

#' @export
#' @rdname defunct
lmtp_ipw <- function(...) {
  lifecycle::deprecate_stop(
    when = "1.5.0",
    what = "lmtp_ipw()",
    details = c("IPW requires the use of correctly pre-specified parametric models for valid statistical inference. Use `lmtp_tmle()` or `lmtp_sdr()` instead.")
  )
}

#' Defunct functions
#'
#' @description
#' `r lifecycle::badge("defunct")`
#'
#' No longer used in the package.
#'
#' @keywords internal
#' @name defunct_funcs
create_node_list <- function(...) {
  lifecycle::deprecate_stop("1.5.0", "create_node_list()")
}
