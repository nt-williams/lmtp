#' LMTP Sequential Doubly Robust Estimator
#'
#' Sequentially doubly robust estimator for the effects of traditional causal effects and
#' modified treatment policies for both point treatment and longitudinal data with binary,
#' continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\] or \[\code{list}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
#'  Or, a list of vectors, the same length as the number of time points of observation.
#'  Vectors should contain column names for the treatment variables at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param outcome \[\code{character}\]\cr
#'  The column name of the outcome variable. In the case of time-to-event
#'  analysis, a vector containing the columns names of intermediate outcome variables and the final
#'  outcome variable ordered by time. Only numeric values are allowed. If the outcome type
#'  is binary, data should be coded as 0 and 1.
#' @param baseline \[\code{character}\]\cr
#'  An optional vector containing the column names of baseline covariates to be
#'  included for adjustment at every time point.
#' @param time_vary \[\code{list}\]\cr
#'  A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param cens \[\code{character}\]\cr
#'  An optional vector of column names of censoring indicators the same
#'  length as the number of time points of observation. If missingness in the outcome is
#'  present or if time-to-event outcome, must be provided.
#' @param shift \[\code{closure}\]\cr
#'  A two argument function that specifies how treatment variables should be shifted.
#'  See examples for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param shifted \[\code{data.frame}\]\cr
#'  An optional data frame, the same as in \code{data}, but modified according
#'  to the treatment policy of interest. If specified, \code{shift} is ignored.
#' @param k \[\code{integer(1)}\]\cr
#'  An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param mtp \[\code{logical(1)}\]\cr
#'  Is the intervention of interest a modified treatment policy?
#'  Default is \code{FALSE}. If treatment variables are continuous this should be \code{TRUE}.
#' @param outcome_type \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., continuous, binomial, survival).
#' @param id \[\code{character(1)}\]\cr
#'  An optional column name containing cluster level identifiers.
#' @param bounds \[\code{numeric(2)}\]\cr
#'  An optional, ordered vector of the bounds for a continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners_outcome \[\code{character}\]\cr A vector of \code{mlr3superlearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"glm"}.
#' @param learners_trt \[\code{character}\]\cr A vector of \code{mlr3superlearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"glm"}.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
#' @param control \[\code{list()}\]\cr
#'  Output of \code{lmtp_control()}.
#'
#' @details
#' ## Should \code{mtp = TRUE}?
#' A modified treatment policy (MTP) is an intervention that depends
#' on the natural value of the exposure (the value that the treatment would have taken under no intervention).
#' This differs from other causal effects,
#' such as the average treatment effect (ATE), where an exposure would be increased (or decreased) deterministically.
#' \bold{If your intervention of interest adds, subtracts, or multiplies the observed treatment values
#' by some amount, use \code{mtp = TRUE}}.
#'
#' @return A list of class \code{lmtp} containing the following components:
#'
#' \item{estimator}{The estimator used, in this case "SDR".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{The estimated standard error of the LMTP effect.}
#' \item{low}{Lower bound of the 95% confidence interval of the LMTP effect.}
#' \item{high}{Upper bound of the 95% confidence interval of the LMTP effect.}
#' \item{eif}{The estimated, un-centered, influence function of the estimate.}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{outcome_reg}{An n x Tau + 1 matrix of outcome regression predictions.
#'  The mean of the first column is used for calculating theta.}
#' \item{density_ratios}{An n x Tau matrix of the estimated, non-cumulative, density ratios.}
#' \item{fits_m}{A list the same length as \code{folds}, containing the fits at each time-point
#'  for each fold for the outcome regression.}
#' \item{fits_r}{A list the same length as \code{folds}, containing the fits at each time-point
#' for each fold of density ratio estimation.}
#' \item{outcome_type}{The outcome variable type.}
#'
#' @example inst/examples/sdr-ex.R
#' @export
lmtp_sdr <- function(data,
                      trt,
                      outcome,
                      baseline = NULL,
                      time_vary = NULL,
                      cens = NULL,
                      shift = NULL,
                      shifted = NULL,
                      k = Inf,
                      mtp = FALSE,
                      outcome_type = c("binomial", "continuous", "survival"),
                      id = NULL,
                      learners_outcome = "glm",
                      learners_trt = "glm",
                      folds = 10,
                      weights = NULL, # NEEDS TO BE CHARACTER NOW
                      control = lmtp_control()) {

  task <- LmtpWideTask$new(
    data = data,
    shifted = make_shifted(data, trt, cens, shift, shifted),
    A = trt,
    Y = outcome,
    W = baseline,
    L = time_vary,
    C = cens,
    id = id,
    weights = weights,
    outcome_type = match.arg(outcome_type),
    mtp = mtp,
    folds = folds
  )

  pb <- progressr::progressor(task$tau*folds*2)

  density_ratios <- crossfit_density_ratio(task, learners_trt, control, pb)
  outcome_regs <- crossfit_sdr(task, density_ratios$ratios, learners_outcome, control, pb)

  theta_dr(task = task,
           m = list(natural = outcome_regs$natural,
                    shifted = outcome_regs$shifted),
           r = density_ratios$ratios,
           fits_m = outcome_regs$fits,
           fits_r = density_ratios$fits,
           shift = deparse(substitute((shift))),
           augmented = FALSE)
}
