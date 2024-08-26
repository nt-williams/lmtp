theta_sub <- function(eta) {
  if (is.null(eta$weights)) {
    theta <- mean(eta$m[, 1])
  }

  if (!is.null(eta$weights)) {
    theta <- weighted.mean(eta$m[, 1], eta$weights)
  }

  if (eta$outcome_type == "continuous") {
    theta <- rescale_y_continuous(theta, eta$bounds)
  }

  out <- list(
    estimator = "substitution",
    theta = theta,
    standard_error = NA_real_,
    low = NA_real_,
    high = NA_real_,
    shift = eta$shift,
    outcome_reg = switch(
      eta$outcome_type,
      continuous = rescale_y_continuous(eta$m, eta$bounds),
      binomial = eta$m
    ),
    fits_m = eta$fits_m,
    outcome_type = eta$outcome_type
  )

  class(out) <- "lmtp"
  out
}

theta_ipw <- function(eta) {
  if (is.null(eta$weights)) {
    theta <- mean(eta$r[, eta$tau]*missing_outcome(eta$y))
  }

  if (!is.null(eta$weights)) {
    theta <- weighted.mean(
      eta$r[, eta$tau]*missing_outcome(eta$y),
      eta$weights
    )
  }

  out <- list(
    estimator = "IPW",
    theta = theta,
    standard_error = NA_real_,
    low = NA_real_,
    high = NA_real_,
    shift = eta$shift,
    density_ratios = eta$r,
    fits_r = eta$fits_r
  )

  class(out) <- "lmtp"
  out
}

eif <- function(r, tau, shifted, natural) {
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, 2:(tau + 1), drop = FALSE] - natural[, 1:tau, drop = FALSE]
  rowSums(compute_weights(r, 1, tau) * m, na.rm = TRUE) + shifted[, 1]
}

theta_dr <- function(eta, augmented = FALSE) {
  inflnce <- eif(r = eta$r,
                 tau = eta$tau,
                 shifted = eta$m$shifted,
                 natural = eta$m$natural)
  theta <- {
    if (augmented)
      weighted.mean(inflnce, eta$weights)
    else
      weighted.mean(eta$m$shifted[, 1], eta$weights)
  }

  if (eta$outcome_type == "continuous") {
    inflnce <- rescale_y_continuous(inflnce, eta$bounds)
    theta <- rescale_y_continuous(theta, eta$bounds)
  }

  clusters <- split(inflnce*eta$weights, eta$id)
  j <- length(clusters)
  se <- sqrt(var(vapply(clusters, function(x) mean(x), 1)) / j)
  ci_low  <- theta - (qnorm(0.975) * se)
  ci_high <- theta + (qnorm(0.975) * se)

  out <- list(
    estimator = eta$estimator,
    theta = theta,
    standard_error = se,
    low = ci_low,
    high = ci_high,
    eif = inflnce,
    id = eta$id,
    shift = eta$shift,
    outcome_reg = switch(
      eta$outcome_type,
      continuous = rescale_y_continuous(eta$m$shifted, eta$bounds),
      binomial = eta$m$shifted
    ),
    density_ratios = eta$r,
    weights = eta$weights,
    fits_m = eta$fits_m,
    fits_r = eta$fits_r,
    outcome_type = eta$outcome_type
  )

  class(out) <- "lmtp"
  out
}
