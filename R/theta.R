
compute_theta <- function(eta, estimator) {

  out <- switch(estimator,
                "ipw" = theta_ipw(r = eta$r, y = eta$y, tau = eta$tau),
                "sub" = theta_sub(m = eta$m, outcome_type = eta$outcome_type, bounds = eta$bounds),
                "tml" = theta_tml_sdr(estimator = "TMLE", m = eta$m, r = eta$r, tau = eta$tau, outcome_type = eta$outcome_type, bounds = eta$bounds, shift = eta$shift),
                "sdr" = theta_tml_sdr(estimator = "SDR", m = eta$m, r = eta$r, tau = eta$tau, outcome_type = eta$outcome_type, bounds = eta$bounds, shift = eta$shift))

  return(out)
}

theta_sub <- function(m, outcome_type, bounds = NULL) {

  # rescale if necessary
  if (outcome_type == "continuous") {
    m <- rescale_y_continuous(m, bounds)
  }

  # calculate estimates
  theta <- mean(m[, 1])

  # returns
  out <- list(estimator = "substitution",
              theta = theta,
              standard_error = NA_real_,
              low = NA_real_,
              high = NA_real_)

  class(out) <- "lmtp"

  return(out)
}

theta_ipw <- function(r, y, tau) {
  # calculate estimates
  theta <- mean(r[, tau]*y, na.rm = T)

  # returns
  out <- list(estimator = "IPW",
              theta = theta,
              standard_error = NA_real_,
              low = NA_real_,
              high = NA_real_)

  class(out) <- "lmtp"

  return(out)
}

eif <- function(r, tau, shifted, natural) {
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, 2:(tau + 1), drop = FALSE] - natural[, 1:tau, drop = FALSE]
  out <- rowSums(r * m, na.rm = TRUE) + shifted[, 1]
  return(out)
}

theta_tml_sdr <- function(estimator, m, r, tau, outcome_type, bounds = NULL,
                          shift) {

  # calculate eif
  inflnce <- eif(r = r, tau = tau, shifted = m$shifted, natural = m$natural)

  # rescale if necessary and calculate estimates
  theta <- mean(m$shifted[, 1])

  if (outcome_type == "continuous") {
    theta <- rescale_y_continuous(theta, bounds)
  }

  n <- nrow(m$natural)
  se <- sd(inflnce, na.rm = TRUE) / sqrt(n)
  ci_low <- theta - (qnorm(0.975) * se)
  ci_high <- theta + (qnorm(0.975) * se)

  # returns
  out <- list(estimator = estimator,
              theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              eif = inflnce,
              shift = shift)

  class(out) <- "lmtp"

  return(out)
}


