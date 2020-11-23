compute_theta <- function(estimator, eta) {
  out <- switch(estimator,
                "ipw" = theta_ipw(eta),
                "sub" = theta_sub(eta),
                "dr"  = theta_dr(eta))

  return(out)
}

theta_sub <- function(eta) {
  i <- Reduce(c, lapply(eta$folds, function(x) x[["validation_set"]]))
  out <- list(estimator = "substitution",
              theta = mean(eta$m[, 1]),
              standard_error = NA_real_,
              low = NA_real_,
              high = NA_real_,
              shift = eta$shift,
              outcome_reg = switch(
                eta$outcome_type,
                continuous = rescale_y_continuous(eta$m, eta$bounds)[order(i), ],
                binomial = eta$m[order(i)]
              ),
              weights_m = eta$weights_m,
              outcome_type = eta$outcome_type)

  class(out) <- "lmtp"
  return(out)
}

theta_ipw <- function(eta) {
  i <- Reduce(c, lapply(eta$folds, function(x) x[["validation_set"]]))
  theta <- mean(eta$r[, eta$tau]*eta$y[i], na.rm = T)

  out <- list(estimator = "IPW",
              theta = theta,
              standard_error = NA_real_,
              low = NA_real_,
              high = NA_real_,
              shift = eta$shift,
              density_ratios = eta$r[order(i), ],
              weights_r = eta$weights_r)

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

theta_dr <- function(eta) {
  i <- Reduce(c, lapply(eta$folds, function(x) x[["validation_set"]]))
  inflnce <- eif(r = eta$r, tau = eta$tau, shifted = eta$m$shifted,
                 natural = eta$m$natural)[order(i)]
  theta <- mean(eta$m$shifted[, 1])

  if (eta$outcome_type == "continuous") {
    inflnce <- rescale_y_continuous(inflnce, eta$bounds)
    theta   <- rescale_y_continuous(theta, eta$bounds)
  }

  clusters <- split(inflnce, eta$id)
  j <- length(clusters)
  se <- sqrt(var(vapply(clusters, function(x) sum(x), 1)) / j)
  ci_low  <- theta - (qnorm(0.975) * se)
  ci_high <- theta + (qnorm(0.975) * se)

  out <- list(estimator = eta$estimator,
              theta = theta,
              standard_error = se,
              low = ci_low,
              high = ci_high,
              eif = inflnce,
              shift = eta$shift,
              outcome_reg = switch(
                eta$outcome_type,
                continuous = rescale_y_continuous(eta$m$shifted, eta$bounds)[order(i), ],
                binomial = eta$m$shifted[order(i), ]
              ),
              density_ratios = eta$r,
              weights_m = eta$weights_m,
              weights_r = eta$weights_r,
              outcome_type = eta$outcome_type)

  class(out) <- "lmtp"
  return(out)
}
