eif <- function(density_ratios, sporadic_weights, shifted, natural, t, tau, l) {
  if (missing(tau)) tau <- ncol(density_ratios)
  if (missing(t)) t <- 1
  if (missing(l)) l <- NULL

  # natural[is.na(natural)] <- -999
  # shifted[is.na(shifted)] <- -999

  # Handle missing sporadic weights
  if (missing(sporadic_weights)) {
    sporadic_weights <- NULL
  }

  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  weights <- compute_weights(density_ratios, sporadic_weights, t, tau, l)
  rowSums(weights * m, na.rm = TRUE) + shifted[, t]
}

compute_weights <- function(density_ratios, sporadic_weights, t, tau, l) {
  # If sporadic_weights is NULL, not in curve algorithm and set to 1
  if (is.null(sporadic_weights)) {
    ipw_sporadic <- 1
  }

  # We only use the sporadic weights if in the first loop of the curve algorithm, l = 1
  else if (l > 1 || is.null(l)) {
    ipw_sporadic <- 1
  }

  else {
    ipw_sporadic <- sporadic_weights[, t:tau, drop = FALSE]
  }

  out <- t(apply(density_ratios[, t:tau, drop = FALSE]*ipw_sporadic, 1, cumprod))
  if (ncol(out) > ncol(density_ratios)) return(t(out))
  out
}
