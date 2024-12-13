eif <- function(r, shifted, natural, t, tau) {
  if (missing(tau)) tau <- ncol(r)
  if (missing(t)) t <- 1
  # natural[is.na(natural)] <- -999
  # shifted[is.na(shifted)] <- -999
  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  rowSums(compute_weights(r, t, tau) * m, na.rm = TRUE) + shifted[, t]
}

compute_weights <- function(r, t, tau) {
  out <- t(apply(r[, t:tau, drop = FALSE], 1, cumprod))
  if (ncol(out) > ncol(r)) return(t(out))
  out
}
