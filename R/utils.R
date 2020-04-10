
.onAttach <- function(libname, pkg) {
  packageStartupMessage(welcome_msg(), check_for_sl3())
}

shift_data <- function(data, A, .f) {

  if (is.null(.f)) {
    return(data)
  }

  out <- data
  sapply(A, function(x) {
    out[, x] <<- .f(data[[x]])
  }, simplify = TRUE)
  return(out)
}

bound <- function(x, p = 1e-5) {
  pmax(pmin(x, 1 - p), p)
}

truncate <- function(x, p = 1e-2) {
  pmin(x, 1 - p)
}

scale_y_continuous <- function(y, bounds) {
  out <- (y - bounds[1]) / (bounds[2] - bounds[1])
  if (is.null(bounds)) {
    out <- y
  }
  return(out)
}

y_bounds <- function(y, outcome_type, bounds = NULL) {
  if (outcome_type == "binomial" || is.null(outcome_type)) {
    out <- NULL
  } else if (is.null(bounds)) {
    out <- c(min(y), max(y))
  } else {
    out <- c(bounds[1], bounds[2])
  }

  return(out)
}

rescale_y_continuous <- function(scaled, bounds) {
  out <- (scaled*(bounds[2] - bounds[1])) + bounds[1]
  return(out)
}

add_scaled_y <- function(data, scaled) {
  data$xyz <- scaled
  return(data)
}

run_ensemble <- function(ensemble, task) {
  ensemble$train(task)
}

predict_sl3 <- function(object, task) {
  out <- object$predict(task)
  return(out)
}

create_censoring_indicators <- function(data, C, tau) {

  # when no censoring return TRUE for all obs
  if (is.null(C)) {
    i <- rep(TRUE, nrow(data))
    j <- rep(TRUE, nrow(data))
    out <- list(i = i, j = j)
    return(out)
  }

  # other wise find censored observations
  i <- data[[C[tau]]] == 1

  if (tau > 1) {
    j <- data[[C[tau - 1]]] == 1
  } else {
    j <- rep(TRUE, nrow(data))
  }

  out <- list(i = i, j = j)
  return(out)
}

