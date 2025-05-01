make_shifted <- function(data, trt, cens, shift, shifted) {
  assert_function(shift, nargs = 2, null.ok = TRUE)

  if (!is.null(shifted)) {
    assert_correctly_shifted(data, shifted, trt, cens)
    return(shifted)
  }

  if (is.null(shifted) && !is.null(shift)) {
    return(shift_data(data, trt, cens, shift))
  }

  if (is.null(shifted) && is.null(shift)) {
    return(shift_data(data, trt, cens, shift))
  }
}

shift_data <- function(data, trt, cens, shift) {
  if (is.null(shift)) {
    return(shift_cens(data, cens))
  }

  is_multivariate <- is.list(trt)
  if (isTRUE(is_multivariate)) {
    return(shift_trt_list(shift_cens(data, cens), trt, shift))
  }

  shift_trt_character(shift_cens(data, cens), trt, shift)
}

shift_cens <- function(data, cens) {
  out <- as.list(data)
  for (ce in cens) {
    out[[ce]] <- 1
  }
  as.data.frame(out, check.names = FALSE)
}

shift_trt_character <- function(data, trt, .f) {
  out <- as.list(data)
  for (a in trt) {
    out[[a]] <- .f(data, a)
  }
  as.data.frame(out, check.names = FALSE)
}

shift_trt_list <- function(data, trt, .f) {
  out <- as.list(data)
  for (a in trt) {
    new <- .f(data, a)
    for (col in a) {
      out[[col]] <- new[[col]]
    }
  }
  as.data.frame(out, check.names = FALSE)
}

#' Turn All Treatment Nodes On
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 1.
#'
#' @param data A dataframe containing the treatment variables.
#' @param trt The name of the current treatment variable.
#'
#' @seealso [lmtp_tmle()], [lmtp_sdr()]
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

#' @seealso [lmtp_tmle()], [lmtp_sdr()]
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

#' IPSI Function Factory
#'
#' A function factory that returns a shift function for increasing or decreasing
#' the probability of exposure when exposure is binary.
#'
#' @param delta \[\code{numeric(1)}\]\cr
#'  A risk ratio between 0 and Inf.
#'
#' @seealso [lmtp_tmle()], [lmtp_sdr()]
#' @return A shift function.
#' @export
#'
#' @examples
#' \donttest{
#' data("iptwExWide", package = "twang")
#' a <- paste0("tx", 1:3)
#' baseline <- c("gender", "age")
#' tv <- list(c("use0"), c("use1"), c("use2"))
#' lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
#'          shift = ipsi(0.5), outcome_type = "continuous", folds = 2)
#' }
ipsi <- function(delta) {
  if (delta > 1) {
    return(ipsi_up(1 / delta))
  }
  ipsi_down(delta)
}

ipsi_up <- function(delta) {
  function(data, trt) {
    eps <- runif(nrow(data), 0, 1)
    ifelse(eps < delta, data[[trt]], 1)
  }
}

ipsi_down <- function(delta) {
  function(data, trt) {
    eps <- runif(nrow(data), 0, 1)
    ifelse(eps < delta, data[[trt]], 0)
  }
}
