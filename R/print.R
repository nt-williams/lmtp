
welcome_msg <- function() {
  cli::cli_text("{.strong lmtp}: Causal Effects Based on Longitudinal Modified Treatment Policies")
  cli::cli_text("{.strong Version}: ", as.character(utils::packageVersion("lmtp")))
}

#' @export
print.lmtp <- function(x, ...) {
  cat("\n")
  cli::cli_text("{.strong LMTP Estimator}: {x$estimator}")
  cli::cli_text(cat("   "), "{.strong Trt. Policy}: ", cli::col_blue(cli::style_italic("{x$shift}")))
  cat("\n")
  cli::cli_text("{.strong Population intervention effect}")
  cli::cli_text(cat("      "), "{.strong Estimate}: {round(x$theta, 4)}")
  cli::cli_text(cat("    "), "{.strong Std. error}: {round(x$standard_error, 4)}")
  cli::cli_text(cat("        "), "{.strong 95% CI}: ({round(x$low, 4)}, {round(x$high, 4)})")
  if (x$estimator %in% c("substitution", "IPW")) no_stderr_warning(x$estimator)
  cat("\n")
}

#' @export
print.lmtp_contrast <- function(x, ...) {
  cat("\n")
  cli::cli_text("{.strong LMTP Contrast}: {x$contrast}")
  cat("\n")
  cli::cli_text(cat("     "), "{.strong Estimate}: {round(x$theta, 4)}")
  cli::cli_text(cat("   "), "{.strong Std. error}: {round(x$standard_error, 4)}")
  cli::cli_text(cat("       "), "{.strong 95% CI}: ({round(x$low, 4)}, {round(x$high, 4)})")
  cli::cli_text(cat("      "), "{.strong p-value}: {format.pval(x$p, eps = 0.00001)}")
  cat("\n")
}
