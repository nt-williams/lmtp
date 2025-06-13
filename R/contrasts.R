#' Perform Contrasts of LMTP Fits
#'
#' Estimates contrasts of multiple LMTP fits compared to either a known reference value
#' or a reference LMTP fit.
#'
#' @param ... One or more objects of class lmtp.
#' @param ref A reference value or another lmtp fit to compare all other fits against.
#' @param type The contrasts of interest. Options are "additive" (the default),
#'  "rr", and "or".
#'
#' @return A list of class \code{lmtp_contrast} containing the following components:
#'
#' \item{type}{The type of contrast performed.}
#' \item{null}{The null hypothesis.}
#' \item{estimates}{A dataframe containing the contrasts estimates, standard errors, and confidence intervals.}
#' @export
#'
#' @example inst/examples/contrasts-ex.R
lmtp_contrast <- function(..., ref, type = c("additive", "rr", "or")) {
  fits <- list(...)
  type <- match.arg(type)

  assert_lmtp_list(fits)
  assert_ref_class(ref)
  assert_contrast_type(type, fits, .var.name = "type")

  if (checkmate::test_number(ref)) {
    cli::cli_alert_danger("Unless {.code ref} is a known constant, standard errors may be anti-conservative!")
  }

  fits <- lapply(fits, function(x) x$estimate)
  if (is.lmtp(ref)) {
    ref <- ref$estimate
  }

  ans <- switch(type,
                "additive" = lapply(fits, function(x) x - ref),
                "rr" = lapply(fits, function(x) x / ref),
                "or" = lapply(fits, function(x) (x / (1 - x)) / (ref / (1 - ref))))

  ans <- do.call("rbind", lapply(ans, function(x) ife::tidy(x)))
  ans$p.value <- pnorm(abs(ans$estimate) / ans$std.error, lower.tail = FALSE) * 2
  ans$shift <- sapply(fits, function(x) x@x)
  ans$ref <- ifelse(inherits(ref, "ife::influence_func_estimate"), ref@x, ref)
  ans <- ans[, c("shift", "ref", "estimate", "std.error", "conf.low", "conf.high", "p.value")]

  structure(list(type = type,
                 null = ifelse(type == "additive", 0, 1),
                 estimates = ans),
            class = "lmtp_contrast")
}
