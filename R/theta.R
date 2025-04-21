theta_dr <- function(task, sequential_regressions, density_ratios, fits_m, fits_r, shift, augmented = FALSE) {
  ic <- eif(density_ratios, shifted = sequential_regressions$shifted, natural = sequential_regressions$natural)

  if (augmented) {
    theta <- weighted.mean(ic, task$weights)
  } else {
    theta <- weighted.mean(sequential_regressions$shifted[, 1], task$weights)
  }

  ic <- task$rescale(ic)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(augmented, "SDR", "TMLE"),
    estimate = ife::ife(theta, ic, task$weights, as.character(task$id)),
    shift = shift,
    outcome_reg = task$rescale(sequential_regressions$shifted),
    density_ratios = density_ratios,
    fits_m = fits_m,
    fits_r = fits_r,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}

theta_curve <- function(task, m, r, sporadic_weights, fits_m, fits_r, shift) {
  ics <- lapply(1:task$tau, function(t) {
    eif(r[, 1:t, drop = FALSE], sporadic_weights[, 1:t, drop = FALSE], m$shifted[[t]], m$natural[[t]])
  })

  thetas <- unlist(lapply(ics, \(x) weighted.mean(x, task$weights)))

  ics <- lapply(ics, \(x) task$rescale(x))
  thetas <- task$rescale(thetas)

  if (task$survival) {
    iso_projection <- isotone::gpava(seq_along(thetas), 1 - thetas)
    thetas <- 1 - iso_projection$x[seq_along(thetas)]
  }

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
