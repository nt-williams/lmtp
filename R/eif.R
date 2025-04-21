eif <- function(density_ratios, sporadic_weights, shifted, natural, t, tau) {
  if (missing(tau)) tau <- ncol(density_ratios)
  if (missing(t)) t <- 1
  # natural[is.na(natural)] <- -999
  # shifted[is.na(shifted)] <- -999
  if (missing(sporadic_weights)) {
    sporadic_weights <- matrix(1, nrow = nrow(density_ratios), ncol = length(t:tau))
  }
  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  rowSums(compute_weights(density_ratios, sporadic_weights, t, tau) * m, na.rm = TRUE) + shifted[, t]
}

compute_weights <- function(density_ratios, sporadic_weights, t, tau) {
  # NEED TO ASK IVAN ABOUT THIS
  gZ <- sporadic_weights[, (t:tau), drop = FALSE]^(c(rep(0, length(t:tau) - 1), 1))
  out <- t(apply(density_ratios[, t:tau, drop = FALSE]*gZ, 1, cumprod))
  if (ncol(out) > ncol(density_ratios)) return(t(out))
  out
}
