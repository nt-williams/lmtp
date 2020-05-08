
.onAttach <- function(libname, pkg) {
  packageStartupMessage(welcome_msg(), check_for_sl3())
}

shift_data <- function(data, A, .f) {
  if (is.null(.f)) {
    return(data)
  }

  out <- as.list(data)
  for (a in A) {
    out[[a]] <- .f(out[[a]])
  }

  return(as.data.frame(out))
}

set_lmtp_options <- function(option, val) {

  if (option == "bound") {
    options(lmtp.bound = val)
  } else if (option == "trunc") {
    option(lmtp.trunc = val)
  }

}

bound <- function(x, p = getOption("lmtp.bound")) {
  pmax(pmin(x, 1 - p), p)
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

transform_sdr <- function(r, tau, max, shifted, natural) {
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, (tau + 2):(max + 1), drop = FALSE] - natural[, (tau + 1):max, drop = FALSE]
  out <- rowSums(r * m, na.rm = TRUE) + shifted[, tau + 1]
  return(out)
}

recombine_ipw <- function(r) {
  out <- list(r = Reduce(rbind, Reduce(rbind, lapply(r, function(x) x[["valid"]]))[, "natural"]),
              sl_weights = lapply(r, function(x) x[["sl_weights"]]))
  return(out)
}

recombine_dens_ratio <- function(r) {
  return(Reduce(rbind, lapply(r, function(x) x[["valid"]])))
}

create_lrnr_matrix <- function(folds, tau, lrnrs) {
  # create empty matrices
  lrnrs <- ifelse(lrnrs == 0, 2, lrnrs)
  out <- lapply(1:folds, function(x) matrix(nrow = tau, ncol = lrnrs))

  # set matrix names
  names(out) <- paste0("fold_", 1:folds)
  for (i in 1:length(out)) {
    rownames(out[[i]]) <- paste0("time_", 1:tau)
    colnames(out[[i]]) <- paste0("Lrnr_", 1:lrnrs)
  }

  return(out)
}

extract_sl_weights <- function(fit) {
  fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE$coefficients
}

count_lrnrs <- function(lrnr) {

  if (is.null(lrnr)) {
    return(2)
  }

  out <- length(lapply(lrnr$params$learners, function(x) x$name))
  if (out == 0) {
    out <- 1
  }
  return(out)
}

pluck_weights <- function(type, x) {
  switch(type,
         "m" = x$sl_weights,
         "r" = lapply(x, function(x) x$sl_weights))
}

is.lmtp <- function(x) {
  class(x) == "lmtp"
}
