#' Longitudinal Targeted Maximum Likelihood Estimator
#'
#' Targeted maximum likelihood estimator for the treatment specific means for
#' discrete exposures with binary, continuous, or time-to-event outcomes.
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\] or \[\code{list}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
#'  Treatment variables should be discrete.
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
#' @param compete \[\code{character}\]\cr
#'  An optional vector of column names of competing risk indicators the same
#'  length as the number of time points of observation. Only used when \code{outcome_type = "survival"}.
#' @param k \[\code{integer(1)}\]\cr
#'  An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param outcome_type \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., continuous, binomial, survival).
#' @param id \[\code{character(1)}\]\cr
#'  An optional column name containing cluster level identifiers.
#' @param bounds \[\code{numeric(2)}\]\cr
#'  An optional, ordered vector of the bounds for a continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners_outcome \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}.
#' @param learners_trt \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param learners_cens \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the censoring mechanism. Default is \code{"SL.glm"}.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
#' @param control \[\code{list()}\]\cr
#'  Output of \code{lmtp_control()}.
#'
#' @details
#'
#' @return A list of class \code{TODO} containing the following components:
#'
#' @example
#' @export
ltmle <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                  cens = NULL, compete = NULL,
                  k = Inf, mtp = TRUE,
                  outcome_type = c("binomial", "continuous", "survival"),
                  id = NULL, bounds = NULL,
                  learners_outcome = "SL.glm",
                  learners_trt = "SL.glm",
                  learners_cens = "SL.glm",
                  folds = 10, weights = NULL,
                  control = lmtp_control()) {

  assert_not_data_table(data)
  variable_names <- c(unlist(trt), outcome, unlist(time_vary), baseline, cens, compete, id)
  assert_subset(variable_names, names(data))
  assert_outcome_types(data, outcome, match.arg(outcome_type))
  assert_numeric(bounds, len = 2, unique = TRUE, sorted = TRUE, finite = TRUE, null.ok = TRUE)
  assert_trt_discrete(data, unlist(trt))

  task <- LmtpTask$new(
    data = data,
    shifted = NULL,
    A = trt,
    Y = outcome,
    L = time_vary,
    W = baseline,
    C = cens,
    D = compete,
    k = k,
    id = id,
    outcome_type = match.arg(outcome_type),
    folds = folds,
    weights = weights,
    bounds = bounds
  )

}
