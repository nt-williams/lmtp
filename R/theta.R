theta_lmtp <- function(task, sequential_regressions, density_ratios, fits_m, fits_r, shift, is_sdr) {
  influence_function <- eif(density_ratios, sequential_regressions$shifted, sequential_regressions$natural)

  if (is_sdr) {
    theta <- fmean(influence_function, w = task$weights)
  } else {
    theta <- fmean(sequential_regressions$shifted[, 1], w = task$weights)
  }

  influence_function <- task$rescale(influence_function)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(is_sdr, "SDR", "TMLE"),
    estimate = ife(theta, influence_function, task$weights, as.character(task$id)),
    shift = shift,
    outcome_reg = task$rescale(sequential_regressions$shifted),
    density_ratios = density_ratios,
    fits_outcome = fits_m,
    fits_treatment = fits_r,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}
