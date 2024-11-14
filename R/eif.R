eif <- function(x, ...) {
  UseMethod("eif")
}

#' @export
eif.matrix <- function(r, shifted, natural) {
  tau <- ncol(r)
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, 2:(tau + 1), drop = FALSE] - natural[, 1:tau, drop = FALSE]
  rowSums(compute_weights(r, 1, tau) * m, na.rm = TRUE) + shifted[, 1]
}

transform_sdr <- function(r, tau, max, shifted, natural) {
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, (tau + 2):(max + 1), drop = FALSE] - natural[, (tau + 1):max, drop = FALSE]
  rowSums(r * m, na.rm = TRUE) + shifted[, tau + 1]
}
