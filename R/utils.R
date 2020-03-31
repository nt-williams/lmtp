
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

rexpit <- function(x) {
  out <- exp(x) / (1 + exp(x))
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

run_ensemble <- function(ensemble) {
  ensemble$stack$train(ensemble$task)
}

predict_sl3_nondensity <- function(object, task) {
  out <- object$predict(task)
  return(out)
}

predict_sl3_density <- function(object, task) {
  out <- object$predict(task)$likelihood
  return(out)
}

theta_sub <- function(m, outcome_type, bounds = NULL, method) {
  if (outcome_type == "continuous") {
    rescaled <- rescale_y_continuous(m, bounds)
    out <- mean(rescaled)
  } else if (outcome_type == "binomial" & method == "glm") {
    out <- mean(rexpit(m))
  } else if (outcome_type == "binomial" & method == "sl") {
    out <- mean(m)
  }
  return(out)
}

theta_ipw <- function(r, y, tau) {
  out <- mean(r[, tau]*y)
  return(out)
}

compute_theta <- function(eta, estimator, outcome_type, bounds = NULL, method = NULL) {

  # TODO: as the rest of the estimators are established need to write their theta methods
  out <- switch(estimator,
                "sub" = theta_sub(m = eta[, 1], outcome_type = outcome_type, bounds = bounds, method = method),
                "ipw" = theta_ipw(r = eta$r, y = eta$y, tau = eta$tau))

  return(out)
}
