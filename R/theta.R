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
