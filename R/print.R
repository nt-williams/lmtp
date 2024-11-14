#' @export
print.lmtp <- function(x, ...) {
  cat("\n")
  cli::cli_text("{.strong LMTP Estimator}: {x$estimator}")
  cli::cli_text(cat("   "), "{.strong Trt. Policy}: ", cli::col_blue(cli::style_italic("{x$shift}")))
  cat("\n")
  cli::cli_text("{.strong Population intervention estimate}")
  print(x$estimate)
  if (x$estimator %in% c("substitution", "IPW")) no_stderr_warning(x$estimator)
  cat("\n")
}

#' @export
print.lmtp_contrast <- function(x, ...) {
  cat("\n")
  cli::cli_text(cat("  "), "{.strong LMTP Contrast}: {x$type}")
  cli::cli_text("{.strong Null hypothesis}: theta == {x$null}")
  cat("\n")
  x$vals$p.value <- format.pval(x$vals$p.value, digits = 3, eps = 0.001)
  print(format(x$vals, digits = 3))
}

#' @export
print.lmtp_survival <- function(x, ...) {
  print(as.data.frame(tidy.lmtp_survival(x)))
}

#' @export
print.lmtp_curve <- function(x, ...) {
  print(as.data.frame(tidy.lmtp_curve(x)))
}
