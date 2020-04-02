
compute_theta <- function(eta, estimator, outcome_type, bounds = NULL, method = NULL) {

  out <- switch(estimator,
                "sub" = theta_sub(m = eta$m, outcome_type = eta$outcome_type, bounds = eta$bounds, method = eta$method),
                "ipw" = theta_ipw(r = eta$r, y = eta$y, tau = eta$tau),
                "tml" = theta_tml_sdr(m = eta$m, outcome_type = eta$outcome_type, bounds = eta$bounds),
                "sdr" = theta_tml_sdr(m = eta$m, outcome_type = eta$outcome_type, bounds = eta$bounds))

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

theta_tml_sdr <- function(m, outcome_type, bounds = NULL) {
  if (outcome_type == "continuous") {
    rescaled <- rescale_y_continuous(m, bounds)
    out <- mean(rescaled)
  } else {
    out <- mean(m)
  }
  return(out)
}


