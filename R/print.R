#' @export
print.lmtp <- function(x, ...) {
  cat("\n")
  cli::cli_text("{.strong LMTP Estimator}: {x$estimator}")
  cli::cli_text(cat("   "), "{.strong Trt. Policy}: ", cli::col_blue(cli::style_italic("{x$shift}")))
  cli::cli_h2("{.emph Population intervention estimate}")
  print(x$estimate)
}

#' @export
print.lmtp_contrast <- function(x, ...) {
  if (nrow(x$estimates) > 1) {
    cli::cli_alert_danger("P-values are not adjusted for multiple comparisons!")
  }
  cat("\n")
  cli::cli_text(cat("  "), "{.strong LMTP Contrast}: {x$type}")
  cli::cli_text("{.strong Null hypothesis}: theta == {x$null}")
  cat("\n")
  x$estimates$p.value <- format.pval(x$estimates$p.value, digits = 3, eps = 0.001)
  print(format(as.data.frame(x$estimates), digits = 3))
}

#' @export
print.lmtp_survival <- function(x, ...) {
  cat("\n")
  cli::cli_text("{.strong LMTP Estimator}: {x[[1]]$estimator}")
  cli::cli_text(cat("   "), "{.strong Trt. Policy}: ", cli::col_blue(cli::style_italic("{x[[1]]$shift}")))
  cli::cli_h2("{.emph Population intervention estimates}")
  print(format(as.data.frame(tidy.lmtp_survival(x)), digits = 3))
}
