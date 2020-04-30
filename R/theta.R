
compute_theta <- function(estimator, eta) {
  out <- switch(estimator,
                "ipw" = theta_ipw(eta),
                "sub" = theta_sub(eta),
                "tml" = theta_tml(eta),
                "sdr" = theta_sdr(eta))

  return(out)
}

theta_sub <- function(eta) {

  # rescale if necessary
  if (eta$outcome_type == "continuous") {
    eta$m <- rescale_y_continuous(eta$m, eta$bounds)
  }

  # calculate estimates
  theta <- mean(eta$m[, 1])

  # returns
  out <- list(estimator = "substitution",
              theta = theta,
              standard_error = NA_real_,
              low = NA_real_,
              high = NA_real_,
              shift = eta$shift)

  class(out) <- "lmtp"

  return(out)
}

theta_ipw <- function(eta) {
  # calculate estimates
  i <- Reduce(c, lapply(eta$folds, function(x) x[["validation_set"]]))
  theta <- mean(eta$r[, eta$tau]*eta$y[i], na.rm = T)

  # returns
  out <- list(estimator = "IPW",
              theta = theta,
              standard_error = NA_real_,
              low = NA_real_,
              high = NA_real_,
              shift = eta$shift)

  class(out) <- "lmtp"

  return(out)
}

eif <- function(r, tau, shifted, natural) {
  natural[is.na(natural)] <- 0
  shifted[is.na(shifted)] <- 0
  m <- shifted[, 2:(tau + 1), drop = FALSE] - natural[, 1:tau, drop = FALSE]
  out <- rowSums(r * m, na.rm = TRUE) + shifted[, 1]
  return(out)
}

theta_tml <- function(eta) {

  i <- Reduce(c, lapply(eta$folds, function(x) x[["validation_set"]]))

  # calculate eif
  inflnce <- eif(r = eta$r, tau = eta$tau,
                 shifted = eta$m$shifted, natural = eta$m$natural)

  # rescale if necessary and calculate estimates
  theta <- mean(eta$m$shifted[, 1])

  if (eta$outcome_type == "continuous") {
    theta <- rescale_y_continuous(theta, eta$bounds)
  }

  n       <- nrow(eta$m$natural)
  se      <- sd(inflnce, na.rm = TRUE) / sqrt(n)
  ci_low  <- theta - (qnorm(0.975) * se)
  ci_high <- theta + (qnorm(0.975) * se)

  # returns
  out <- list(estimator = "TMLE",
              theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              eif = inflnce[order(i)],
              shift = eta$shift)

  class(out) <- "lmtp"

  return(out)
}

theta_sdr <- function(eta) {

  i <- Reduce(c, lapply(eta$folds, function(x) x[["validation_set"]]))

  # calculate eif
  inflnce <- eif(r = eta$r, tau = eta$tau,
                 shifted = eta$m$shifted, natural = eta$m$natural)

  # rescale if necessary and calculate estimates
  theta <- mean(eta$m$shifted[, 1])

  if (eta$outcome_type == "continuous") {
    inflnce <- rescale_y_continuous(inflnce, eta$bounds)
    theta   <- rescale_y_continuous(theta, eta$bounds)
  }

  n       <- nrow(eta$m$natural)
  se      <- sd(inflnce, na.rm = TRUE) / sqrt(n)
  ci_low  <- theta - (qnorm(0.975) * se)
  ci_high <- theta + (qnorm(0.975) * se)

  # returns
  out <- list(estimator = "SDR",
              theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              eif = inflnce[order(i)],
              shift = eta$shift)

  class(out) <- "lmtp"

  return(out)
}
