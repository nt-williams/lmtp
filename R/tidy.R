#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a(n) lmtp object
#'
#' @param x A `lmtp` object produced by a call to [lmtp::lmtp_tmle()], [lmtp::lmtp_sdr()],
#' [lmtp::lmtp_survival()].
#' @param ... Unused, included for generic consistency only.
#'
#' @examples
#' \donttest{
#' a <- c("A1", "A2")
#' nodes <- list(c("L1"), c("L2"))
#' cens <- c("C1", "C2")
#' y <- "Y"
#' fit <- lmtp_tmle(sim_cens, a, y, time_vary = nodes, cens = cens, shift = NULL, folds = 2)
#' tidy(fit)
#' }
#'
#' @export
tidy.lmtp <- function(x, ...) ife::tidy(x$estimate)

#' Tidy a(n) lmtp_survival object
#'
#' @param x A `lmtp_survival` object produced by a call to [lmtp::lmtp_survival()].
#' @param ... Unused, included for generic consistency only.
#'
#' @example inst/examples/lmtp_survival-ex.R
#'
#' @export
tidy.lmtp_survival <- function(x, ...) {
  out <- do.call("rbind", lapply(x, tidy))
  out$time <- seq_along(x)
  out[, c(ncol(out), 1:ncol(out) - 1)]
}
