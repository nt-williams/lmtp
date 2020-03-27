
shift_data <- function(data, A, shift) {
  out <- data
  sapply(A, function(x) {
    out[, x] <<- data[[x]] + shift
  }, simplify = TRUE)
  return(out)
}

create_m <- function(n, t, Y) {
  out <- matrix(nrow = n, ncol = t)
  out[, t] <- Y
  return(out)
}

rexpit <- function(x) {
  out <- exp(x) / (1 + exp(x))
  return(out)
}

scale_y_continuous <- function(x, mi = NULL, ma = NULL) {

  if (is.null(mi)) {
    mi <- min(x)
  }

  if (is.null(ma)) {
    ma <- max(x)
  }

  out <- (x - mi) / (ma - mi)
  return(out)
}

run_ensemble <- function(ensemble) {
  ensemble$stack$train(ensemble$task)
}
