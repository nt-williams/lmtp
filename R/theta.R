theta_dr <- function(task, sequential_regressions, density_ratios, fits_m, fits_r, shift, is_sdr) {
  influence_function <- eif(density_ratios, sequential_regressions$shifted, sequential_regressions$natural)

  if (is_sdr) {
    theta <- weighted.mean(influence_function, task$weights)
  } else {
    theta <- weighted.mean(sequential_regressions$shifted[, 1], task$weights)
  }

  influence_function <- task$rescale(influence_function)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(is_sdr, "SDR", "TMLE"),
    estimate = ife::ife(theta, influence_function, task$weights, as.character(task$id)),
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

theta_ltmle <- function(task, estimates, propensity_scores, levels) {
  theta <- sapply(estimates, function(x) fmean(x$predictions[, 1], w = task$weights))

  # Rescale estimates
  theta <- sapply(theta, function(x) task$rescale(x))
  influence_functions <- lapply(estimates, function(x) task$rescale(x$uncentered_eif))

  # Create 'ife' objects
  ifes <- lapply(seq_along(theta), function(i) {
    ife::ife(theta[i], influence_functions[[i]], task$weights, as.character(task$id))
  })
  names(estimates) <- levels

  out <- list(
    estimator = "TMLE",
    estimates = ifes,
    outcome_reg = lapply(estimates, function(x) task$rescale(x$predictions)),
    propensity_scores = propensity_scores$propensity_score,
    prob_observed = propensity_scores$prob_observed,
    fits_outcome = lapply(estimates, function(x) x$learner_outcome_summary),
    fits_treatment = propensity_scores$learner_treatment_summary,
    fits_censoring = propensity_scores$learner_cens_summary,
    outcome_type = task$outcome_type
  )

  class(out) <- "ltmle"
  out
}
