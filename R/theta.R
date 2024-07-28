theta_sub <- function(eta) {
  cumulative_indicator <- as.logical(apply(eta$conditional, 1, prod))

  if (is.null(eta$weights)) {
    theta <- mean(eta$m[cumulative_indicator, 1])
  }

  if (!is.null(eta$weights)) {
    theta <- weighted.mean(eta$m[cumulative_indicator, 1], eta$weights[cumulative_indicator])
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
  cumulative_indicator <- as.logical(apply(eta$conditional, 1, prod))

  if (is.null(eta$weights)) {
    theta <- mean(eta$r[, eta$tau]*missing_outcome(eta$y)) / mean(cumulative_indicator)
  }

  if (!is.null(eta$weights)) {
    theta <- weighted.mean(
      eta$r[, eta$tau]*missing_outcome(eta$y),
      eta$weights
    ) / weighted.mean(cumulative_indicator, eta$weights)
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

eif <- function(r, cumulated, tau, shifted, natural, conditional, G) {
  cumulative_indicator <- as.logical(apply(conditional, 1, prod))
  future_indicator <- t(apply(conditional[, ncol(conditional):1], 1, cumprod))[,ncol(conditional):1]

  natural[is.na(natural)] <- 0
  shifted[is.na(shifted)] <- 0
  m <- shifted[, 2:(tau + 1), drop = FALSE] - natural[, 1:tau, drop = FALSE]

  if (cumulated == TRUE) {
    weights <- r
  } else {
    weights <- compute_weights(r, 1, tau)
  }
  theta <- mean(shifted[cumulative_indicator, 1])
  #1 / mean(cumulative_indicator) * (rowSums(weights * m, na.rm = TRUE) + cumulative_indicator * (shifted[, 1] - theta))
  theta + rowSums(weights * future_indicator[,2:(tau + 1)] / G * m, na.rm = TRUE) / mean(cumulative_indicator)
}

theta_dr <- function(eta, augmented = FALSE) {
  cumulative_indicator <- as.logical(apply(eta$conditional, 1, prod))

  inflnce <- eif(r = eta$r,
                 cumulated = eta$cumulated,
                 tau = eta$tau,
                 shifted = eta$m$shifted,
                 natural = eta$m$natural,
                 conditional = eta$conditional,
                 G = eta$G)

  theta <- {
    if (augmented)
      weighted.mean(eta$m$shifted[cumulative_indicator, 1], eta$weights[cumulative_indicator]) # why isn't this using 'inflnce'?
    else
      weighted.mean(eta$m$shifted[cumulative_indicator, 1], eta$weights[cumulative_indicator])
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
