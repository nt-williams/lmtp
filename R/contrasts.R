
#' Estimate intervention contrast effects
#'
#' @param x An object of class lmtp.
#' @param y An object of class lmtp.
#' @param contrast The contrast to perform between \code{x} and \code{y}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
lmtp_contrast <- function(x, y,
                          contrast = c("additive", "relative risk", "odds ratio")) {

  out <- switch(match.arg(contrast),
                "additive" = contrast_additive(x = x, y = y),
                "relative risk" = contrast_relative_risk(x = x, y = y),
                "odds ratio" = contrast_odds_ratio(x = x, y = y))

  class(out) <- "lmtp_contrast"

  print(x)
  print(y)
  return(out)
}

contrast_additive <- function(x, y) {
  theta <- x$theta - y$theta
  eif <- x$eif - y$eif
  se <- sd(eif) / sqrt(length(eif))
  ci_low <- theta - qnorm(0.975) * se
  ci_high <- theta + qnorm(0.975) * se
  p <- pnorm(abs(theta) / se, lower.tail = FALSE) * 2

  out <- list(theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              p = p,
              contrast = "Additive intervention effect",
              eif = eif)

  return(out)
}

contrast_relative_risk <- function(x, y) {
  theta <- x$theta / y$theta
  log_eif <- (x$eif / x$theta) - (y$eif / y$theta)
  se <- sd(log_eif) / sqrt(length(log_eif))
  ci_low <- exp(log(theta) - qnorm(0.975) * se)
  ci_high <- exp(log(theta) + qnorm(0.975) * se)
  p <- pnorm(abs(log(theta)) / se, lower.tail = FALSE) * 2

  out <- list(theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              p = p,
              contrast = "Relative risk",
              eif = eif)

  return(out)
}

contrast_odds_ratio <- function(x, y) {
  theta <- (x$theta / (1 - x$theta)) / (y$theta / (1 - y$theta))
  log_eif <- (x$eif / (x$theta * (1 - x$theta))) - (y$eif / (y$theta * (1 - y$theta)))
  se <- sd(log_eif) / sqrt(length(log_eif))
  ci_low <- exp(log(theta) - qnorm(0.975) * se)
  ci_high <- exp(log(theta) + qnorm(0.975) * se)
  p <- pnorm(abs(log(theta)) / se, lower.tail = FALSE) * 2

  out <- list(theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              p = p,
              contrast = "Odds ratio",
              eif = eif)

  return(out)
}
