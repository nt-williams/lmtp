
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
  cat("\n")
  cli::cli_text("Recommended package, {.pkg sl3}, not detected.")
  cli::cli_text("{.pkg sl3} can be installed with: {.code remotes::install_github('tlverse/sl3')}")
  cat("\n")
}

check_pb <- function(pb, t, status) {
  if (isFALSE(pb) | t == 1) {
    pb <- NULL
  } else {
    pb <- initiate_progress_bar(status, tau = t)
  }
  return(pb)
}

check_sd <- function(x, learner_stack) {
  if (sd(x) > .Machine$double.eps) {
    out <- learner_stack
  } else {
    out <- sl3::Lrnr_mean
  }

  return(out)
}


