eif <- function(density_ratios, sporadic_weights, shifted, natural, t, tau, l) {
  if (missing(tau)) tau <- ncol(density_ratios)
  if (missing(t)) t <- 1
  # If missing l, then we are not in the curve algorithm and sporadic weights don't apply
  if (missing(l)) l <- 1

  # If missing sporadic weights, then we assume all weights are 1
  if (missing(sporadic_weights)) {
    sporadic_weights <- matrix(1, nrow = nrow(density_ratios), ncol = length(t:tau))
  }

  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  weights <- compute_weights(density_ratios, sporadic_weights, t, tau, l)
  rowSums(weights * m, na.rm = TRUE) + shifted[, t]
}

compute_weights <- function(density_ratios, sporadic_weights, t, tau, l) {
  # We only use the sporadic weights if in the first loop of the curve algorithm
  ipw_sporadic <- sporadic_weights[, (t:tau), drop = FALSE]^(as.numeric(l == 1))
  out <- t(apply(density_ratios[, t:tau, drop = FALSE]*ipw_sporadic, 1, cumprod))
  if (ncol(out) > ncol(density_ratios)) return(t(out))
  out
}
