
.onLoad <- function(libname, pkg) {
  packageStartupMessage(welcome_msg(), check_for_sl3())
}

welcome_msg <- function() {
  cat("\n")
  cli::cli_text("{.strong lmtp}: Causal Effects Based on Longitudinal Modified Treatment Policies")
  cli::cli_text("{.strong Version}: ", as.character(packageVersion("lmtp")))
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
