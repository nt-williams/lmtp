#' LMTP Targeted Maximum Likelihood Estimator
#'
#' Targeted maximum likelihood estimator for the effects of traditional causal effects and
#' modified treatment policies for both point treatment and longitudinal data with binary,
#' continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
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
#' @param learners_outcome \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param learners_trt \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
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
#'  Return full SuperLearner fits? Default is \code{FALSE}, return only SuperLearner weights.
#' @param ... Extra arguments. Exists for backwards compatibility.
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
#' \item{estimator}{The estimator used, in this case "TMLE".}
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
#' @example inst/examples/tmle-ex.R
#' @export
lmtp_tmle <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                      cens = NULL, shift = NULL, shifted = NULL, k = Inf,
                      mtp = FALSE, outcome_type = c("binomial", "continuous", "survival"),
                      # intervention_type = c("static", "dynamic", "mtp"),
                      id = NULL, bounds = NULL,
                      learners_outcome = "SL.glm",
                      learners_trt = "SL.glm",
                      folds = 10, weights = NULL, .bound = 1e-5, .trim = 0.999,
                      .learners_outcome_folds = 10, .learners_trt_folds = 10,
                      .return_full_fits = FALSE, ...) {

  assertNotDataTable(data)
  checkmate::assertCharacter(outcome, len = if (match.arg(outcome_type) != "survival") 1,
                             min.len = if (match.arg(outcome_type) == "survival") 2)
  checkmate::assertCharacter(baseline, null.ok = TRUE)

  tau <- determine_tau(outcome, trt)

  assertTrtCharacter(trt, tau)
  checkmate::assertCharacter(cens, len = tau, null.ok = !checkmate::anyMissing(data[, outcome, drop = FALSE]))
  checkmate::assertList(time_vary, types = c("NULL", "character"), len = tau, null.ok = TRUE)
  checkmate::assertCharacter(id, len = 1, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  assertLmtpData(data, trt, outcome, baseline, time_vary, cens, id)
  assertOutcomeTypes(data, outcome, match.arg(outcome_type))
  assertReservedNames(data)
  checkmate::assertFunction(shift, nargs = 2, null.ok = TRUE)
  assertShiftedData(shifted, data, c(outcome, baseline, unlist(time_vary), id), cens)
  checkmate::assertNumeric(bounds, len = 2, finite = TRUE, any.missing = FALSE, sorted = TRUE, null.ok = TRUE)
  checkmate::assertNumeric(weights, len = nrow(data), finite = TRUE, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumber(k, lower = 0, upper = Inf)
  checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)
  checkmate::assertNumber(.learners_outcome_folds, null.ok = TRUE)
  checkmate::assertNumber(.learners_trt_folds, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  checkmate::assertNumber(.bound)
  checkmate::assertNumber(.trim, upper = 1)
  checkmate::assertLogical(.return_full_fits, len = 1)

  extras <- list(...)
  if ("intervention_type" %in% names(extras)) {
    mtp <- extras$intervention_type == "mtp"
    warning("The `intervention_type` argument of `lmtp_tmle()` is deprecated as of lmtp 1.3.1",
            call. = FALSE)
  }

  Task <- lmtp_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = weights,
    bounds = bounds,
    bound = .bound
  )

  pb <- progressr::progressor(Task$tau*folds*2)

  ratios <- cf_r(Task, learners_trt, mtp, .learners_trt_folds, .trim, .return_full_fits, pb)
  estims <- cf_tmle(Task, "tmp_lmtp_scaled_outcome", ratios$ratios, learners_outcome, .learners_outcome_folds, .return_full_fits, pb)

  theta_dr(
    list(
      estimator = "TMLE",
      m = list(natural = estims$natural, shifted = estims$shifted),
      r = ratios$ratios,
      tau = Task$tau,
      folds = Task$folds,
      id = Task$id,
      outcome_type = Task$outcome_type,
      bounds = Task$bounds,
      weights = Task$weights,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      fits_m = estims$fits,
      fits_r = ratios$fits,
      outcome_type = Task$outcome_type
    ),
    FALSE
  )
}

#' LMTP Sequential Doubly Robust Estimator
#'
#' Sequentially doubly robust estimator for the effects of traditional causal effects and
#' modified treatment policies for both point treatment and longitudinal data with binary,
#' continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
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
#' @param learners_outcome \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param learners_trt \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
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
#'  Return full SuperLearner fits? Default is \code{FALSE}, return only SuperLearner weights.
#' @param ... Extra arguments. Exists for backwards compatibility.
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
lmtp_sdr <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                     cens = NULL, shift = NULL, shifted = NULL, k = Inf,
                     mtp = FALSE,
                     # intervention_type = c("static", "dynamic", "mtp"),
                     outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL,
                     learners_outcome = "SL.glm",
                     learners_trt = "SL.glm",
                     folds = 10, weights = NULL, .bound = 1e-5, .trim = 0.999,
                     .learners_outcome_folds = 10, .learners_trt_folds = 10,
                     .return_full_fits = FALSE, ...) {

  assertNotDataTable(data)
  checkmate::assertCharacter(outcome, len = if (match.arg(outcome_type) != "survival") 1,
                             min.len = if (match.arg(outcome_type) == "survival") 2)
  checkmate::assertCharacter(baseline, null.ok = TRUE)

  tau <- determine_tau(outcome, trt)

  assertTrtCharacter(trt, tau)
  checkmate::assertCharacter(cens, len = tau, null.ok = !checkmate::anyMissing(data[, outcome, drop = FALSE]))
  checkmate::assertList(time_vary, types = c("NULL", "character"), len = tau, null.ok = TRUE)
  checkmate::assertCharacter(id, len = 1, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  assertLmtpData(data, trt, outcome, baseline, time_vary, cens, id)
  assertOutcomeTypes(data, outcome, match.arg(outcome_type))
  assertReservedNames(data)
  checkmate::assertFunction(shift, nargs = 2, null.ok = TRUE)
  assertShiftedData(shifted, data, c(outcome, baseline, unlist(time_vary), id), cens)
  checkmate::assertNumeric(bounds, len = 2, finite = TRUE, any.missing = FALSE, sorted = TRUE, null.ok = TRUE)
  checkmate::assertNumeric(weights, len = nrow(data), finite = TRUE, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumber(k, lower = 0, upper = Inf)
  checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)
  checkmate::assertNumber(.learners_outcome_folds, null.ok = TRUE)
  checkmate::assertNumber(.learners_trt_folds, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  checkmate::assertNumber(.bound)
  checkmate::assertNumber(.trim, upper = 1)
  checkmate::assertLogical(.return_full_fits, len = 1)

  Task <- lmtp_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = weights,
    bounds = bounds,
    bound = .bound
  )

  extras <- list(...)
  if ("intervention_type" %in% names(extras)) {
    mtp <- extras$intervention_type == "mtp"
    warning("The `intervention_type` argument of `lmtp_sdr()` is deprecated as of lmtp 1.3.1",
            call. = FALSE)
  }

  pb <- progressr::progressor(Task$tau*folds*2)

  ratios <- cf_r(Task, learners_trt, mtp, .learners_trt_folds, .trim, .return_full_fits, pb)
  estims <- cf_sdr(Task, "tmp_lmtp_scaled_outcome", ratios$ratios, learners_outcome, .learners_outcome_folds, .return_full_fits, pb)

  theta_dr(
    list(
      estimator = "SDR",
      m = list(natural = estims$natural, shifted = estims$shifted),
      r = ratios$ratios,
      tau = Task$tau,
      folds = Task$folds,
      id = Task$id,
      outcome_type = Task$outcome_type,
      bounds = Task$bounds,
      weights = Task$weights,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      fits_m = estims$fits,
      fits_r = ratios$fits,
      outcome_type = Task$outcome_type
    ),
    TRUE
  )
}

#' LMTP Substitution Estimator
#'
#' G-computation estimator for the effects of traditional causal effects and
#' modified treatment policies for both point treatment and longitudinal data with binary,
#' continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
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
#' @param outcome_type \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., continuous, binomial, survival).
#' @param id \[\code{character(1)}\]\cr
#'  An optional column name containing cluster level identifiers.
#' @param bounds \[\code{numeric(2)}\]\cr
#'  An optional, ordered vector of the bounds for a continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
#' @param .bound \[\code{numeric(1)}\]\cr
#'  Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .learners_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners}.
#' @param .return_full_fits \[\code{logical(1)}\]\cr
#'  Return full SuperLearner fits? Default is \code{FALSE}, return only SuperLearner weights.
#'
#' @return A list of class \code{lmtp} containing the following components:
#'
#' \item{estimator}{The estimator used, in this case "substitution".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{NA}
#' \item{low}{NA}
#' \item{high}{NA}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{outcome_reg}{An n x Tau + 1 matrix of outcome regression predictions.
#'  The mean of the first column is used for calculating theta.}
#' \item{fits_m}{A list the same length as \code{folds}, containing the fits at each time-point
#' for each fold for the outcome regression.}
#' \item{outcome_type}{The outcome variable type.}
#' @export
#'
#' @example inst/examples/sub-ex.R
lmtp_sub <- function(data, trt, outcome, baseline = NULL, time_vary = NULL, cens = NULL,
                     shift = NULL, shifted = NULL, k = Inf,
                     outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL, learners = "SL.glm",
                     folds = 10, weights = NULL, .bound = 1e-5, .learners_folds = 10,
                     .return_full_fits = FALSE) {

  assertNotDataTable(data)
  checkmate::assertCharacter(outcome, len = if (match.arg(outcome_type) != "survival") 1,
                             min.len = if (match.arg(outcome_type) == "survival") 2)
  checkmate::assertCharacter(baseline, null.ok = TRUE)

  tau <- determine_tau(outcome, trt)

  assertTrtCharacter(trt, tau)
  checkmate::assertCharacter(cens, len = tau, null.ok = !checkmate::anyMissing(data[, outcome, drop = FALSE]))
  checkmate::assertList(time_vary, types = c("NULL", "character"), len = tau, null.ok = TRUE)
  checkmate::assertCharacter(id, len = 1, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  assertLmtpData(data, trt, outcome, baseline, time_vary, cens, id)
  assertOutcomeTypes(data, outcome, match.arg(outcome_type))
  assertReservedNames(data)
  checkmate::assertFunction(shift, nargs = 2, null.ok = TRUE)
  assertShiftedData(shifted, data, c(outcome, baseline, unlist(time_vary), id), cens)
  checkmate::assertNumeric(bounds, len = 2, finite = TRUE, any.missing = FALSE, sorted = TRUE, null.ok = TRUE)
  checkmate::assertNumeric(weights, len = nrow(data), finite = TRUE, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumber(k, lower = 0, upper = Inf)
  checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)
  checkmate::assertNumber(.learners_folds, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  checkmate::assertNumber(.bound)
  checkmate::assertLogical(.return_full_fits, len = 1)

  Task <- lmtp_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = weights,
    bounds = bounds,
    bound = .bound
  )

  pb <- progressr::progressor(Task$tau*folds)

  estims <- cf_sub(Task, "tmp_lmtp_scaled_outcome", learners, .learners_folds, .return_full_fits, pb)

  theta_sub(
    eta = list(
      m = estims$m,
      outcome_type = Task$outcome_type,
      bounds = Task$bounds,
      folds = Task$folds,
      weights = Task$weights,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      fits_m = estims$fits,
      outcome_type = Task$outcome_type
    )
  )
}

#' LMTP IPW Estimator
#'
#' Inverse probability of treatment weighting estimator for the effects of traditional causal
#' effects and modified treatment policies for both point treatment and longitudinal data
#' with binary, continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
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
#' @param learners \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
#' @param .bound \[\code{numeric(1)}\]\cr
#'  Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .trim \[\code{numeric(1)}\]\cr
#'  Determines the amount the density ratios should be trimmed.
#'  The default is 0.999, trimming the density ratios greater than the 0.999 percentile
#'  to the 0.999 percentile. A value of 1 indicates no trimming.
#' @param .learners_folds \[\code{integer(1)}\]\cr
#'  The number of cross-validation folds for \code{learners}.
#' @param .return_full_fits \[\code{logical(1)}\]\cr
#'  Return full SuperLearner fits? Default is \code{FALSE}, return only SuperLearner weights.
#' @param ... Extra arguments. Exists for backwards compatibility.
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
#' \item{estimator}{The estimator used, in this case "IPW".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{NA}
#' \item{low}{NA}
#' \item{high}{NA}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{density_ratios}{An n x Tau matrix of the estimated density ratios.}
#' \item{fits_r}{A list the same length as \code{folds}, containing the fits at each time-point
#' for each fold of density ratio estimation.}
#' @export
#'
#' @example inst/examples/ipw-ex.R
lmtp_ipw <- function(data, trt, outcome, baseline = NULL, time_vary = NULL, cens = NULL,
                     shift = NULL, shifted = NULL, mtp = FALSE,
                     # intervention_type = c("static", "dynamic", "mtp"),
                     k = Inf, id = NULL,
                     outcome_type = c("binomial", "continuous", "survival"),
                     learners = "SL.glm",
                     folds = 10, weights = NULL,
                     .bound = 1e-5, .trim = 0.999, .learners_folds = 10,
                     .return_full_fits = FALSE, ...) {

  assertNotDataTable(data)
  checkmate::assertCharacter(outcome, len = if (match.arg(outcome_type) != "survival") 1,
                             min.len = if (match.arg(outcome_type) == "survival") 2)
  checkmate::assertCharacter(baseline, null.ok = TRUE)

  tau <- determine_tau(outcome, trt)

  assertTrtCharacter(trt, tau)
  checkmate::assertCharacter(cens, len = tau, null.ok = !checkmate::anyMissing(data[, outcome, drop = FALSE]))
  checkmate::assertList(time_vary, types = c("NULL", "character"), len = tau, null.ok = TRUE)
  checkmate::assertCharacter(id, len = 1, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  assertLmtpData(data, trt, outcome, baseline, time_vary, cens, id)
  assertOutcomeTypes(data, outcome, match.arg(outcome_type))
  assertReservedNames(data)
  checkmate::assertFunction(shift, nargs = 2, null.ok = TRUE)
  assertShiftedData(shifted, data, c(outcome, baseline, unlist(time_vary), id), cens)
  checkmate::assertNumeric(weights, len = nrow(data), finite = TRUE, any.missing = FALSE, null.ok = TRUE)
  checkmate::assertNumber(k, lower = 0, upper = Inf)
  checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)
  checkmate::assertNumber(.learners_folds, null.ok = TRUE)
  checkmate::assertSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(data))
  checkmate::assertNumber(.bound)
  checkmate::assertNumber(.trim, upper = 1)
  checkmate::assertLogical(.return_full_fits, len = 1)

  Task <- lmtp_Task$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = weights,
    bounds = NULL,
    bound = .bound
  )

  pb <- progressr::progressor(Task$tau*folds)

  extras <- list(...)
  if ("intervention_type" %in% names(extras)) {
    mtp <- extras$intervention_type == "mtp"
    warning("The `intervention_type` argument of `lmtp_ipw()` is deprecated as of lmtp 1.3.1",
            call. = FALSE)
  }

  ratios <- cf_r(Task, learners, mtp, .learners_folds, .trim, .return_full_fits, pb)

  theta_ipw(
    eta = list(
      r = matrix(
        t(apply(ratios$ratios, 1, cumprod)),
        nrow = nrow(ratios$ratios),
        ncol = ncol(ratios$ratios)
      ),
      y = if (Task$survival) {
        convert_to_surv(data[[final_outcome(outcome)]])
      } else {
        data[[final_outcome(outcome)]]
      },
      folds = Task$folds,
      weights = Task$weights,
      tau = Task$tau,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      fits_r = ratios$fits
    )
  )
}
