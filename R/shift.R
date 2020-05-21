
shift_data <- function(data, trt, cens, shift) {
  if (is.null(shift)) {
    return(shift_cens(data, cens))
  }
  return(shift_trt(shift_cens(data, cens), trt, shift))
}

shift_cens <- function(data, cens) {
  out <- as.list(data)
  for (ce in cens) {
    out[[ce]] <- 1
  }
  return(as.data.frame(out))
}

shift_trt <- function(data, trt, .f) {
  for (a in trt) {
    data[[a]] <- .f(data, a)
  }
  return(data)
}

#' Turn All Treatment Nodes On
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 1.
#'
#' @param data A dataframe containg the treatment variables.
#' @param trt The name of the current treatment variable.
#'
#' @seealso [lmtp_tmle()], [lmtp_sdr()], [lmtp_sub()], [lmtp_ipw()]
#' @return
#' @export
#'
#' @examples
#' # TODO
static_binary_on <- function(data, trt) {
  rep(1, length(data[[trt]]))
}

#' Turn All Treatment Nodes Off
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 0.
#'
#' @param data A dataframe containg the treatment variables.
#' @param trt The name of the current treatment variable.

#' @seealso [lmtp_tmle()], [lmtp_sdr()], [lmtp_sub()], [lmtp_ipw()]
#' @return
#' @export
#'
#' @examples
#' #TODO
static_binary_off <- function(data, trt) {
  rep(0, length(data[[trt]]))
}