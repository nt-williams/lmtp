
prepare_mbased <- function(data, A, Y, shift, outcome_type, bounds) {

  check_scaled_conflict(names(data))

  n           <- nrow(data)
  t           <- length(A)
  m           <- cbind(matrix(nrow = n, ncol = t), data[, Y])
  shifted     <- shift_data(data, A, shift)
  scaled      <- scale_y_continuous(data[, Y], outcome_type, bounds)
  shifted$xyz <- data$xyz <- scaled$scaled

  out <- list(n = n,
              tau = t,
              m = m,
              data = data,
              shifted_data = shifted,
              outcome_type = outcome_type,
              scale_meta = scaled)

  return(out)
}

prepare_rbased <- function(data, A, Y, shift) {
  n <- nrow(data)
  t <- length(A)

  out <- list(n = n,
              tau = t)

  return(out)
}
