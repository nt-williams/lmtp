eif <- function(x, ...) {
  UseMethod("eif")
}

#' @export
eif.matrix <- function(r, shifted, natural, t, tau) {
  if (missing(tau)) tau <- ncol(r)
  if (missing(t)) t <- 1
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  rowSums(compute_weights(r, t, tau) * m, na.rm = TRUE) + shifted[, t]
}

transform_sdr <- function(r, tau, max, shifted, natural) {
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, (tau + 2):(max + 1), drop = FALSE] - natural[, (tau + 1):max, drop = FALSE]
  rowSums(r * m, na.rm = TRUE) + shifted[, tau + 1]
}
