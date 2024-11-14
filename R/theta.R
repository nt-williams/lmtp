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

theta_dr <- function(task, m, r, fits_m, fits_r, shift, augmented = FALSE) {
  ic <- eif(r, m$shifted, m$natural)

  if (augmented) {
    theta <- weighted.mean(ic, task$weights)
  } else {
    theta <- weighted.mean(m$shifted[, 1], task$weights)
  }

  ic <- task$rescale(ic)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(augmented, "SDR", "TMLE"),
    estimate = ife::ife(theta, ic, task$weights, as.character(task$id)),
    shift = shift,
    outcome_reg = task$rescale(m$shifted),
    density_ratios = r,
    fits_m = fits_m,
    fits_r = fits_r,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}

theta_curve <- function(task, m, r, fits_m, fits_r, shift) {
  ics <- lapply(1:task$tau, function(t) {
    eif(r[, 1:t, drop = FALSE], m$shifted[[t]], m$natural[[t]])
  })

  thetas <- unlist(lapply(ics, \(x) weighted.mean(x, task$weights)))

  ics <- lapply(ics, \(x) task$rescale(x))
  thetas <- task$rescale(thetas)

  iso_projection <- isotone::gpava(seq_along(thetas), 1 - thetas)
  thetas <- 1 - iso_projection$x[seq_along(thetas)]

  out <- list(
    estimator = "SDR curve",
    estimates = lapply(
      seq_along(ics),
      \(x) ife::ife(thetas[[x]], ics[[x]], task$weights, as.character(task$id))
    ),
    outcome_reg = lapply(m$shifted, \(x) task$rescale(x)),
    density_ratios = r,
    fits_m = fits_m,
    fits_r = fits_r,
    shift = shift,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp_curve"
  out
}

# TODO: NEED TO SAVE THE SEED FOR THE REPLICATES AND THE BOOTED ESTIMATES FOR ESTIMATNG CONTRASTS
theta_boot <- function(eta) {
  theta <- weighted.mean(eta$m[, 1], eta$weights)

  if (eta$outcome_type == "continuous") {
    theta <- rescale_y_continuous(theta, eta$bounds)
  }

  # TODO: NEED TO FIGURE OUT HOW THIS WOULD WORK WITH CLUSTERING
  se <- sqrt(var(eta$boots))
  ci_low  <- theta - (qnorm(0.975) * se)
  ci_high <- theta + (qnorm(0.975) * se)

  out <- list(
    estimator = eta$estimator,
    theta = theta,
    standard_error = se,
    low = ci_low,
    high = ci_high,
    boots = eta$boots,
    id = eta$id,
    shift = eta$shift,
    outcome_reg = switch(
      eta$outcome_type,
      continuous = rescale_y_continuous(eta$m, eta$bounds),
      binomial = eta$m
    ),
    density_ratios = eta$r,
    fits_m = eta$fits_m,
    fits_r = eta$fits_r,
    outcome_type = eta$outcome_type,
    seed = eta$seed
  )

  class(out) <- "lmtp"
  out
}
