eif <- function(density_ratios, sporadic_weights, shifted, natural, t, tau) {
  if (missing(tau)) tau <- ncol(density_ratios)
  if (missing(t)) t <- 1

  # natural[is.na(natural)] <- -999
  # shifted[is.na(shifted)] <- -999

  # Handle missing sporadic weights
  if (missing(sporadic_weights)) {
    sporadic_weights <- NULL
  }

  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  weights <- compute_weights(density_ratios, sporadic_weights, t, tau)
  rowSums(weights * m, na.rm = TRUE) + shifted[, t]
}

compute_weights <- function(density_ratios, sporadic_weights, t, tau) {
  # Sporadic weights only effect the final density ratio
  if (!is.null(sporadic_weights)) {
    density_ratios[, tau] <- density_ratios[, tau] * sporadic_weights[, tau]
  }

  out <- t(apply(density_ratios[, t:tau, drop = FALSE], 1, cumprod))
  if (ncol(out) > ncol(density_ratios)) return(t(out))
  out
}
