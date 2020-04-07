
.onAttach <- function(libname, pkg) {
  packageStartupMessage(welcome_msg(), check_for_sl3())
}

shift_data <- function(data, A, .f) {
  out <- data
  sapply(A, function(x) {
    out[, x] <<- .f(data[[x]])
  }, simplify = TRUE)
  return(out)
}

create_m <- function(n, t, Y) {
  out <- matrix(nrow = n, ncol = t)
  out[, t] <- Y
  return(out)
}

bound <- function(x, p = 1e-5) {
  pmax(pmin(x, 1 - p), p)
}

truncate <- function(x, p = 1e-2) {
  pmin(x, 1 - p)
}

scale_y_continuous <- function(x, outcome_type, bounds = NULL) {
  if (outcome_type == "binomial") {
    out <- list(scaled = x,
                bounds = NULL)
    return(out)
  }

  if (is.null(bounds)) {
    mi <- min(x)
    ma <- max(x)
  } else {
    mi <- bounds[1]
    ma <- bounds[2]
  }

  scaled <- (x - mi) / (ma - mi)
  out <- list(scaled = scaled,
              bounds = c(mi, ma))
  return(out)
}

rescale_y_continuous <- function(scaled, bounds) {
  mi <- bounds[1]
  ma <- bounds[2]

  out <- (scaled*(ma - mi)) + mi
  return(out)
}

run_ensemble <- function(ensemble, task) {
  ensemble$train(task)
}

predict_sl3 <- function(object, task) {
  out <- object$predict(task)
  return(out)
}

create_censoring_indicators <- function(data, C, tau) {
  i <- data[[C[tau]]] == 1

  if (tau > 1) {
    j <- data[[C[tau - 1]]] == 1
  } else {
    j <- rep(TRUE, nrow(data))
  }

  out <- list(i = i,
              j = j)

  return(out)
}
