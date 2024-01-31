#' Set LMTP Estimation Parameters
#'
#' @param .trim \[\code{numeric(1)}\]\cr
#'  Determines the amount the density ratios should be trimmed.
#'  The default is 0.999, trimming the density ratios greater than the 0.999 percentile
#'  to the 0.999 percentile. A value of 1 indicates no trimming.
#' @param .learners_outcome_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners_outcome}.
#' @param .learners_trt_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners_trt}.
#' @param .return_full_fits \[\code{logical(1)}\]\cr
#'  Return full \code{mlr3superlearner} fits? Default is \code{FALSE}.
#'
#' @return A list of parameters controlling the estimation procedure.
#' @export
#'
#' @examples
#' lmtp_control(.trim = 0.975)
lmtp_control <- function(...) {
  change <- list(...)
  control <- list(.trim = 0.999,
                  .learners_outcome_folds = NULL,
                  .learners_trt_folds = NULL,
                  .return_full_fits = FALSE)
  if (length(change) == 0) return(control)
  change <- change[names(change) %in% names(control)]
  control[names(change)] <- change
  control
}
