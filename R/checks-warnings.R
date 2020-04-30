
check_for_sl3 <- function(test = FALSE) {
  tryCatch(
    if (isTRUE(test)) {
      stop()
    } else {
      has_sl3 <- "sl3" %in% rownames(utils::installed.packages())
      if (isFALSE(has_sl3)) stop()
      else on.exit()
    }, error = function(e) {
      no_sl3()
    }
  )
}

no_stderr_warning <- function(estimator) {
  cat("\n")
  cli::cli_alert_warning("Standard errors aren't provided for the {estimator} estimator.")
}

no_sl3 <- function() {
  cli::cli_text("Recommended package, {.pkg sl3}, not detected.")
  cli::cli_text("{.pkg sl3} can be installed with: {.code remotes::install_github('tlverse/sl3')}")
}

check_sd <- function(x, learner_stack) {
  if (sd(x) > .Machine$double.eps) {
    out <- learner_stack
  } else {
    out <- sl3::Lrnr_mean
  }

  return(out)
}

check_censoring <- function(data, training, validation, C, Y, tau) {

  if (any(is.na(data[[Y]])) & is.null(C)) {
    stop("Missing outcomes detected and censoring nodes not indicated.", call. = FALSE)
  } else if (!is.null(C)) {
    check <- TRUE
  } else if (is.null(C) | !any(is.na(data[[Y]]))) {
    check <- FALSE
  }

  ct <- matrix(nrow = nrow(training), ncol = tau)
  cv <- matrix(nrow = nrow(validation), ncol = tau)
  if (isFALSE(check)) {
    for (t in 1:tau) {
      ct[, t] <- rep(1, nrow(training))
      cv[, t] <- rep(1, nrow(validation))
    }
  }

  out <- list(train = ct,
              valid = cv)

  return(out)
}

fix_censoring_ind <- function(data, C = NULL, tau) {

  if (!is.null(C)) {
    out <- as.list(data)
    for (t in 1:tau) {
      out[[C[[t]]]] <- ifelse(is.na(out[[C[[t]]]]), 0, out[[C[[t]]]])
    }
    return(as.data.frame(out))
  }

  return(data)
}

check_scaled_conflict <- function(data) {
  nn <- names(data)
  check <- "xyz" %in% nn

  tryCatch(
    if (check) {
      stop()
    } else {
      on.exit()
    },
    error = function(e) {
      stop(
        "A variable named `xyz` was detected in your data. This variable name is reserved for an interal `lmtp` process. Please rename this variable.",
        call. = FALSE
      )
    }
  )
}

check_extreme_ratio <- function(ratio) {
  return(apply(ratio, 2, function(x) pmin(x, quantile(x, 0.999))))
}

check_variation <- function(data, outcome, learners) {
  if (sd(data[, outcome]) < .Machine$double.eps) {
    learners <- sl3::make_learner(sl3::Lrnr_mean)
  }
  return(learners)
}
