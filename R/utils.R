
.onAttach <- function(libname, pkg) {
  packageStartupMessage(welcome_msg(), check_for_sl3())
}

shift_data <- function(data, A, C, .f) {

  out <- as.list(data)

  if (is.null(.f)) { # only set C = 1
    for (ce in C) {
      out[[ce]] <- 1
    }
  } else {
    for (a in A) { # shift A
      out[[a]] <- .f(out[[a]])
    }

    for (ce in C) { # and set C = 1
      out[[ce]] <- 1
    }
  }

  return(as.data.frame(out))
}

set_lmtp_options <- function(option, val) {
  if (option == "bound") {
    options(lmtp.bound = val)
  } else if (option == "trt") {
    options(lmtp.trt.length = val)
  } else {
    stop("Unknown lmtp option.", call. = F)
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
    out <- c(min(y, na.rm = T), max(y, na.rm = T))
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

create_determ_indicators <- function(data, determ, tau) {
  if (is.null(determ)) {
    return(rep(FALSE, nrow(data)))
  } else {
    return(data[[determ[tau]]] == 1 & !is.na(data[[determ[tau]]]))
  }
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

hold_lrnr_weights <- function(folds) {
  lapply(1:folds, function(x) list())
}

extract_sl_weights <- function(fit) {
  as.data.frame(fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE$fits)
}

pluck_weights <- function(type, x) {
  switch(type,
         "m" = x$sl_weights,
         "r" = lapply(x, function(x) x$sl_weights))
}

is.lmtp <- function(x) {
  class(x) == "lmtp"
}

sw <- function(x) {
  suppressWarnings(x)
}

final_outcome <- function(outcomes) {
  outcomes[length(outcomes)]
}

#' Time To Event Last Outcome Carried Forward
#'
#' A helper function to prepare survival data for use with LMTP estimators
#' by imputing outcome nodes using last outcome carried forward when an observation
#' experiences the event before the end-of-follow-up.
#'
#' @param data The dataset to modify.
#' @param outcomes A vector of outcome nodes ordered by time.
#'
#' @return A modified dataset with future outcome nodes set to 1 if an observation
#'  experienced an event at any previous time point.
#' @export
#' @examples
#' event_locf(sim_point_surv, c(paste0("Y.", 0:6)))
event_locf <- function(data, outcomes) {
  tau <- length(outcomes)
  for (i in outcomes[1:(tau - 1)]) {
    modify <- setdiff(outcomes[match(i, outcomes):tau], i) # find all outcomes at later time
    for (j in 1:nrow(data)) {
      if (data[j, i] == 1 & !is.na(data[j, i])) data[j, modify] <- 1 # if the event is observed at the current time, set all later events to 1
    }
  }
  return(data)
}

