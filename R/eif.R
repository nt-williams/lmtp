eif <- function(density_ratios, shifted, natural, time, time_horizon) {
  if (missing(time_horizon)) time_horizon <- ncol(density_ratios)
  if (missing(time)) time <- 1
  residuals <- shifted[, (time + 1):(time_horizon + 1), drop = FALSE] - natural[, time:time_horizon, drop = FALSE]
  rowSums(compute_weights(density_ratios, time, time_horizon) * residuals, na.rm = TRUE) + shifted[, time]
}

compute_weights <- function(density_ratios, time, time_horizon) {
  out <- t(apply(density_ratios[, time:time_horizon, drop = FALSE], 1, cumprod))
  if (ncol(out) > ncol(density_ratios)) return(t(out))
  out
}

delay_weights <- function(task, propensity, levels, time, time_horizon) {
  # Compare the observed treatment sequence to the current sequence
  D <- apply(task$natural[, task$vars$A[time:time_horizon], drop = FALSE],
             1, function(x) x == levels, simplify = FALSE)
  D <- do.call("rbind", D)
  D <- apply(D, 1, prod)

  # Get the propensity score ratio for each time
  ratios <- sapply(seq_along(time:time_horizon), function(i) {
    t <- (time:time_horizon)[i]
    trt <- task$natural[[task$vars$A[t]]]
    gs <- propensity[[t]][, as.character(levels[i])]
    gA <- this_propensity(trt, propensity[[t]])
    gs / gA
  })

  # ASK IVAN: this just returns a bunch of 0's and 1's
  D * apply(ratios, 1, prod)
}

update_htlmtp_eif <- function(eif, augmented, A, Y, id, predictions,
                              propensity_scores, prob_observed, this_A, this_C, time) {
  # Construct component Riesz representer
  riesz <- delay_riesz_rep(augmented, A, propensity_scores, prob_observed, this_A, this_C, 1, time)
  # Multiply Riesz representer by residual
  component <- riesz * (augmented[[Y]] - predictions)
  # Add to the given EIF
  eif + fsum(component, augmented[[id]])
}
