theta_lmtp <- function(task, sequential_regressions, density_ratios, fits_m, fits_r, shift, is_sdr) {
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

theta_htlmtp <- function(task, estimator, propensity_scores, shift, is_sdr) {
  if (is_sdr) {
    theta <- weighted.mean(estimator$uncentered_eif, task$weights)
  } else {
    theta <- weighted.mean(estimator$predictions, task$weights)
  }

  influence_function <- task$rescale(estimator$uncentered_eif)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(is_sdr, "History of Treatment - SDR", "History of Treatment - TMLE"),
    estimate = ife::ife(theta, influence_function, task$weights, as.character(task$id)),
    shift = shift,
    outcome_reg = task$rescale(estimator$predictions),
    propensity_scores = propensity_scores$propensity_score,
    prob_observed = propensity_scores$prob_observed,
    fits_outcome = estimator$learner_outcome_summary,
    fits_treatment = propensity_scores$learner_treatment_summary,
    fits_censoring = propensity_scores$learner_cens_summary,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}
