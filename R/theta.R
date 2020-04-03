
compute_theta <- function(eta, estimator) {

  out <- switch(estimator,
                "sub" = theta_sub(m = eta$m, outcome_type = eta$outcome_type, bounds = eta$bounds, method = eta$method),
                "ipw" = theta_ipw(r = eta$r, y = eta$y, tau = eta$tau),
                "tml" = theta_tml_sdr(m = eta$m, r = eta$r, tau = eta$tau, outcome_type = eta$outcome_type, bounds = eta$bounds),
                "sdr" = theta_tml_sdr(m = eta$m, r = eta$r, tau = eta$tau, outcome_type = eta$outcome_type, bounds = eta$bounds))

  return(out)
}

eif <- function(r, tau, shifted, natural) {
  m <- shifted[, 2:(tau + 1), drop = FALSE] - natural[, 1:tau, drop = FALSE]
  out <- rowSums(r * m) + shifted[, 1]
  return(out)
}

theta_sub <- function(m, outcome_type, bounds = NULL, method) {
  if (outcome_type == "continuous") {
    rescaled <- rescale_y_continuous(m, bounds)
    out <- mean(rescaled)
  } else if (outcome_type == "binomial" & method == "glm") {
    out <- mean(plogis(m))
  } else if (outcome_type == "binomial" & method == "sl") {
    out <- mean(m)
  }
  return(out)
}

theta_ipw <- function(r, y, tau) {
  out <- mean(r[, tau]*y)
  return(out)
}

theta_tml_sdr <- function(m, r, tau, outcome_type, bounds = NULL) {

  # calculate eif
  inflnce <- eif(r = r, tau = tau, shifted = m$shifted, natural = m$natural)

  # calculate estimates
  theta <- mean(m$shifted[, 1])
  se <- sd(inflnce) / sqrt(length(inflnce))

  # rescale if necessary
  if (outcome_type == "continuous") {
    inflnce <- rescale_y_continuous(inflnce, bounds)
    theta <- rescale_y_continuous(theta, bounds)
    se <- rescale_y_continuous(se, bounds)
  }

  # returns
  out <- list(theta = theta,
              standard_error = se,
              eif = inflnce)

  return(out)
}


