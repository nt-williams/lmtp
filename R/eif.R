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
