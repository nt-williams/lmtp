
check_for_sl3 <- function(test = FALSE) {
  tryCatch(
    if (isTRUE(test)) {
      stop()
    } else {
      has_sl3 <- "sl3" %in% rownames(installed.packages())
      if (isFALSE(has_sl3)) stop()
      else on.exit()
    }, error = function(e) {
      no_sl3()
    }
  )
}

no_stderr_warning <- function(estimator) {
  cli::cli_alert_warning("Standard errors aren't provided for the {estimator} estimator.")
  cat("\n")
}

no_sl3 <- function() {
  cat("\n")
  cli::cli_text("Recommended package, {.pkg sl3}, not detected.")
  cli::cli_text("{.pkg sl3} can be installed with: {.code remotes::install_github('tlverse/sl3')}")
  cat("\n")
}


