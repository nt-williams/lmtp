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
#' \item{vals}{A dataframe containing the contrasts estimates, standard errors, and confidence intervals.}
#' \item{eifs}{Un-centered estimated influence functions for contrasts estimated.}
#' @export
#'
#' @example inst/examples/contrasts-ex.R
lmtp_contrast <- function(..., ref, type = c("additive", "rr", "or")) {
  fits <- list(...)

  check_lmtp_type(fits, ref)
  check_outcome_type(fits, ref, match.arg(type))

  switch(check_ref_type(ref, match.arg(type)),
         "additive" = contrast_additive(fits = fits, ref = ref),
         "rr"       = contrast_rr(fits = fits, ref = ref),
         "or"       = contrast_or(fits = fits, ref = ref))
}

contrast_additive <- function(fits, ref) {
  res <- lapply(fits, function(x) contrast_additive_single(x, ref))
  vals <- Reduce(rbind, lapply(res, function(x) x[["vals"]]))
  eifs <- Reduce(cbind, lapply(res, function(x) x[["eif"]]))

  out <- list(
    type = "additive",
    null = 0,
    vals = vals,
    eifs = eifs
  )
  class(out) <- "lmtp_contrast"
  return(out)
}

contrast_additive_single <- function(fit, ref) {
  if (is.lmtp(ref)) {
    theta <- fit$theta - ref$theta
    eif <- fit$eif - ref$eif
  }

  if (isFALSE(is.lmtp(ref))) {
    theta <- fit$theta - ref
    eif <- fit$eif
  }

  std.error <- sd(eif) / sqrt(length(eif))
  conf.low <- theta - qnorm(0.975) * std.error
  conf.high <- theta + qnorm(0.975) * std.error
  p.value <- pnorm(abs(theta) / std.error, lower.tail = FALSE) * 2

  list(
    vals = data.frame(
      theta = theta,
      shift = fit$theta,
      ref = ifelse(is.lmtp(ref), ref$theta, ref),
      std.error = std.error,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value
    ),
    eif = eif
  )
}

contrast_rr <- function(fits, ref) {
  res <- lapply(fits, function(x) contrast_rr_single(x, ref))
  vals <- Reduce(rbind, lapply(res, function(x) x[["vals"]]))
  eifs <- Reduce(cbind, lapply(res, function(x) x[["eif"]]))

  out <- list(
    type = "relative risk",
    null = 1,
    vals = vals,
    eifs = eifs
  )

  class(out) <- "lmtp_contrast"
  return(out)
}

contrast_rr_single <- function(fit, ref) {
  theta <- fit$theta / ref$theta
  log_eif <- (fit$eif / fit$theta) - (ref$eif / ref$theta)
  std.error <- sd(log_eif) / sqrt(length(log_eif))
  conf.low <- exp(log(theta) - qnorm(0.975) * std.error)
  conf.high <- exp(log(theta) + qnorm(0.975) * std.error)
  p.value <- pnorm(abs(log(theta)) / std.error, lower.tail = FALSE) * 2

  list(
    vals = data.frame(
      theta = theta,
      shift = fit$theta,
      ref = ref$theta,
      std.error = std.error,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value
    ),
    eif = log_eif
  )
}

contrast_or <- function(fits, ref) {
  res <- lapply(fits, function(x) contrast_or_single(x, ref))
  vals <- Reduce(rbind, lapply(res, function(x) x[["vals"]]))
  eifs <- Reduce(cbind, lapply(res, function(x) x[["eif"]]))

  out <- list(
    type = "odds ratio",
    null = 1,
    vals = vals,
    eifs = eifs
  )
  class(out) <- "lmtp_contrast"
  out
}

contrast_or_single <- function(fit, ref) {
  theta <- (fit$theta / (1 - fit$theta)) / (ref$theta / (1 - ref$theta))
  log_eif <- (fit$eif / (fit$theta * (1 - fit$theta))) - (ref$eif / (ref$theta * (1 - ref$theta)))
  std.error <- sd(log_eif) / sqrt(length(log_eif))
  conf.low  <- exp(log(theta) - qnorm(0.975) * std.error)
  conf.high <- exp(log(theta) + qnorm(0.975) * std.error)
  p.value <- pnorm(abs(log(theta)) / std.error, lower.tail = FALSE) * 2

  list(
    vals = data.frame(
      theta = theta,
      shift = fit$theta,
      ref = ref$theta,
      std.error = std.error,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value
    ),
    eif = log_eif
  )
}

