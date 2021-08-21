shift_data <- function(data, trt, cens, shift) {
  if (is.null(shift)) {
    return(shift_cens(data, cens))
  }
  shift_trt(shift_cens(data, cens), trt, shift)
}

shift_cens <- function(data, cens) {
  out <- as.list(data)
  for (ce in cens) {
    out[[ce]] <- 1
  }
  as.data.frame(out, check.names = FALSE)
}

shift_trt <- function(data, trt, .f) {
  for (a in trt) {
    data[[a]] <- .f(data, a)
  }
  data
}

#' Turn All Treatment Nodes On
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 1.
#'
#' @param data A dataframe containing the treatment variables.
#' @param trt The name of the current treatment variable.
#'
#' @seealso [lmtp_tmle()], [lmtp_sdr()], [lmtp_sub()], [lmtp_ipw()]
#' @return A dataframe with all treatment nodes set to 1.
#' @export
#'
#' @examples
#' \donttest{
#' data("iptwExWide", package = "twang")
#' a <- paste0("tx", 1:3)
#' baseline <- c("gender", "age")
#' tv <- list(c("use0"), c("use1"), c("use2"))
#' lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
#'          shift = static_binary_on, outcome_type = "continuous", folds = 2)
#' }
static_binary_on <- function(data, trt) {
  rep(1, length(data[[trt]]))
}

#' Turn All Treatment Nodes Off
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 0.
#'
#' @param data A dataframe containing the treatment variables.
#' @param trt The name of the current treatment variable.

#' @seealso [lmtp_tmle()], [lmtp_sdr()], [lmtp_sub()], [lmtp_ipw()]
#' @return A dataframe with all treatment nodes set to 0.
#' @export
#'
#' @examples
#' \donttest{
#' data("iptwExWide", package = "twang")
#' a <- paste0("tx", 1:3)
#' baseline <- c("gender", "age")
#' tv <- list(c("use0"), c("use1"), c("use2"))
#' lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
#'          shift = static_binary_off, outcome_type = "continuous", folds = 2)
#' }
static_binary_off <- function(data, trt) {
  rep(0, length(data[[trt]]))
}
