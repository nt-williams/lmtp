theta_lmtp <- function(task, estimates, density_ratios, shift, is_sdr) {
  if (is_sdr) {
    theta <- fmean(estimates$uncentered_eif, w = task$weights)
  } else {
    theta <- fmean(estimates$predictions[, 1], w = task$weights)
  }

  influence_function <- task$rescale(estimates$uncentered_eif)
  theta <- task$rescale(theta)

  out <- list(
    estimator = ifelse(is_sdr, "SDR", "TMLE"),
    estimate = ife(theta, influence_function, task$weights, as.character(task$id)),
    shift = shift,
    outcome_reg = task$rescale(estimates$predictions),
    density_ratios = density_ratios$density_ratios,
    fits_outcome = estimates$learner_outcome_summary,
    fits_treatment = density_ratios$learner_treatment_summary,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}

theta_ltmle <- function(task, estimates, propensity_scores, levels, trt_balance, cens_balance) {
  theta <- sapply(estimates, function(x) fmean(x$predictions[, 1], w = task$weights))

  # Rescale estimates
  theta <- sapply(theta, function(x) task$rescale(x))
  influence_functions <- lapply(estimates, function(x) task$rescale(x$uncentered_eif))

  # Create 'ife' objects
  ifes <- lapply(seq_along(theta), function(i) {
    ife::ife(theta[i], influence_functions[[i]], task$weights, as.character(task$id))
  })
  names(ifes) <- levels

  out <- list(
    estimator = "TMLE",
    estimates = ifes,
    outcome_reg = lapply(estimates, function(x) task$rescale(x$predictions)),
    propensity_scores = propensity_scores$propensity_score,
    prob_observed = propensity_scores$prob_observed,
    balance = list(
      treatment = setNames(trt_balance, paste0("time ", seq_len(length(task$vars$A)))),
      censoring = if (!is.null(cens_balance)) setNames(cens_balance, paste0("time ", seq_len(task$time_horizon))) else NULL
    ),
    fits_outcome = lapply(estimates, function(x) x$learner_outcome_summary),
    fits_treatment = propensity_scores$learner_treatment_summary,
    fits_censoring = propensity_scores$learner_cens_summary,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp_ltmle"
  out
}
