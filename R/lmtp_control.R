#' Set LMTP Estimation Parameters
#'
#' @param .bound \[\code{numeric(1)}\]\cr
#'  Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .trim \[\code{numeric(1)}\]\cr
#'  Determines the amount the density ratios should be trimmed.
#'  The default is 0.999, trimming the density ratios greater than the 0.999 percentile
#'  to the 0.999 percentile. A value of 1 indicates no trimming.
#' @param .learners_outcome_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners_outcome}.
#' @param .learners_trt_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners_trt}.
#' @param .return_full_fits \[\code{logical(1)}\]\cr
#'  Return full 'SuperLearner' fits? Default is \code{FALSE}, return only 'SuperLearner' coefficients.
#' @param .discrete \[\code{logical(1)}\]\cr
#'  Use discrete or ensemble super learner?
#' @param .info \[\code{logical(1)}\]\cr
#'  Print super learner fitting info to the console?
#'
#' @return A list of parameters controlling the estimation procedure.
#' @export
#'
#' @examples
#' lmtp_control(.trim = 0.975)
lmtp_control <- function(.bound = 1e5,
                         .trim = 0.999,
                         .learners_outcome_folds = 10,
                         .learners_trt_folds = 10,
                         .return_full_fits = FALSE,
                         .discrete = TRUE,
                         .info = FALSE) {

  assert_number(.learners_outcome_folds, null.ok = TRUE)
  assert_number(.learners_trt_folds, null.ok = TRUE)
  assert_number(.bound)
  assert_number(.trim, upper = 1)
  assert_logical(.return_full_fits, len = 1)
  assert_logical(.discrete, len = 1)
  assert_logical(.info, len = 1)

  list(.bound = .bound,
       .trim = .trim,
       .learners_outcome_folds = .learners_outcome_folds,
       .learners_trt_folds = .learners_trt_folds,
       .return_full_fits = .return_full_fits,
       .discrete = .discrete,
       .info = .info)
}
