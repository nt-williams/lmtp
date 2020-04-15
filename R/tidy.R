
#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a(n) lmtp object
#'
#' @param x A `lmtp` object produced by a call to [lmtp::lmtp_tmle()], [lmtp::lmtp_sdr()],
#' [lmtp::lmtp_sub()], or [lmtp::lmtp_ipw()].
#' @param ... Unused, included for generic consistency only.
#'
#' @examples # TODO
#'
#' @export
tidy.lmtp <- function(x, ...) {
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::tibble(estimator = x$estimator,
                   estimate = x$theta,
                   std.error = x$standard_error,
                   conf.low = x$low,
                   conf.high = x$high)
  } else {
    as.data.frame(estimator = x$estimator,
                  estimate = x$theta,
                  std.error = x$standard_error,
                  conf.low = x$low,
                  conf.high = x$high)
  }
}
