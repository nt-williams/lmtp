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
#' @param .learners_conditional_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners_conditional}.
#' @param .return_full_fits \[\code{logical(1)}\]\cr
#'  Return full SuperLearner fits? Default is \code{FALSE}, return only SuperLearner weights.
#' @param .epochs \[\code{integer(1)}\]\cr
#' The number of epochs to train the neural network.
#' @param .learning_rate \[\code{numeric(1)}\]\cr
#'  The learning rate for the neural network.
#' @param .batch_size \[\code{integer(1)}\]\cr
#'  The batch size for the neural network.
#' @param .device \[\code{character(1)}\]\cr
#'  The device to train the neural network on. Default is \code{"cpu"}.
#'
#' @return A list of parameters controlling the estimation procedure.
#' @export
#'
#' @examples
#' lmtp_control(.trim = 0.975)
lmtp_control <- function(.bound = 1e5,
                         .trim = 0.999,
                         .learners_outcome_folds = NULL,
                         .learners_trt_folds = NULL,
                         .learners_conditional_folds = NULL,
                         .return_full_fits = FALSE,
                         .epochs = 500L,
                         .learning_rate = 0.1,
                         .batch_size = 64,
                         .patience = 10,
                         .weight_decay = 0.001) {
  list(.bound = .bound,
       .trim = .trim,
       .learners_outcome_folds = .learners_outcome_folds,
       .learners_trt_folds = .learners_trt_folds,
       .learners_conditional_folds = .learners_conditional_folds,
       .return_full_fits = .return_full_fits,
       .epochs = .epochs,
       .learning_rate = .learning_rate,
       .batch_size = .batch_size,
       .patience = .patience,
       .weight_decay = .weight_decay)
}
