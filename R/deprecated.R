#' LMTP Substitution Estimator
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been deprecated. Please use [lmtp_tmle()] or [lmtp_sdr()] instead.
#'
#' @keywords internal
#'
#' @param ... Ignored
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' \dontrun{
#' # This function is deprecated
#' # Use lmtp_tmle() or lmtp_sdr() instead
#' }
#'
#' @export
lmtp_sub <- function(...) {
  lifecycle::deprecate_stop(
    when = "2.0.0",
    what = "lmtp_sub()",
    details = c("It requires the use of correctly pre-specified parametric models for valid statistical inference. Use `lmtp_tmle()` or `lmtp_sdr()`.",
                "Visit <https://beyondtheate.com/info_estimators.html> for more information")
  )
}

#' LMTP IPW Estimator
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been deprecated. Please use [lmtp_tmle()] or [lmtp_sdr()] instead.
#'
#' @keywords internal
#'
#' @param ... Ignored
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' \dontrun{
#' # This function is deprecated
#' # Use lmtp_tmle() or lmtp_sdr() instead
#' }
#'
#' @export
lmtp_ipw <- function(...) {
  lifecycle::deprecate_stop(
    when = "2.0.0",
    what = "lmtp_sub()",
    details = c("It requires the use of correctly pre-specified parametric models for valid statistical inference. Use `lmtp_tmle()` or `lmtp_sdr()`.",
                "Visit <https://beyondtheate.com/info_estimators.html> for more information")
  )
}
