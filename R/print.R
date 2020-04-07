
welcome_msg <- function() {
  cli::cli_text("{.strong lmtp}: Causal Effects Based on Longitudinal Modified Treatment Policies")
  cli::cli_text("{.strong Version}: ", as.character(utils::packageVersion("lmtp")))
}

#' @export
print.lmtp <- function(x, ...) {
  cat("\n")
  cli::cli_text("{.strong LMTP Estimator}: {x$estimator}")
  cat("\n")
  cli::cli_text("{.strong Population intervention effect}")
  cli::cli_text(cat("   "), "{.strong Estimate}: {round(x$theta, 4)}")
  cli::cli_text(cat(" "), "{.strong Std. error}: {round(x$standard_error, 4)}")
  cli::cli_text(cat("     "), "{.strong 95% CI}: ({round(x$low, 4)}, {round(x$high, 4)})")
  if (x$estimator %in% c("substitution", "IPW")) no_stderr_warning(x$estimator)
  cat("\n")
}

# progress bar
initiate_progress_bar <- function(section, tau) {

  form <- paste(" ", cli::col_white(section), cli::col_white("[:bar] :percent :elapsed"))
  pb <- progress::progress_bar$new(
    format = form,
    total = tau,
    clear = FALSE,
    show_after = 0.05
  )

  return(pb)
}

progress_progress_bar <- function(pb) {
  if (is.null(pb)) {
    on.exit()
  } else {
    pb$tick()
  }
}
