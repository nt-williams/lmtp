#' LMTP Targeted Maximum Likelihood Estimator
#'
#' Targeted maximum likelihood estimator for the effects of traditional causal effects and
#' modified treatment policies for both point treatment and longitudinal data with binary,
#' continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data A data frame in wide format containing all necessary variables
#'  for the estimation problem.
#' @param trt A vector containing the column names of treatment variables ordered by time.
#' @param outcome The column name of the outcome variable. In the case of time-to-event
#'  analysis, a vector containing the columns names of intermediate outcome variables and the final
#'  outcome variable ordered by time. Only numeric values are allowed. If the outcome type
#'  is binary, data should be coded as 0 and 1.
#' @param baseline An optional vector containing the column names of baseline covariates to be
#'  included for adjustment at every time point.
#' @param time_vary A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as the number of time points of observation. If missingness in the outcome is
#'  present or if time-to-event outcome, must be provided.
#' @param shift A two argument function that specifies how treatment variables should be shifted.
#'  See examples for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param shifted An optional data frame, the same as in \code{data}, but modified according
#'  to the treatment policy of interest. If specified, \code{shift} is ignored.
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param intervention_type The intervention type, should be one of \code{"static"},
#'   \code{"dynamic"}, \code{"mtp"}.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial, survival).
#' @param id An optional column name containing cluster level identifiers.
#' @param bounds An optional vector of the bounds for continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners_outcome A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param learners_trt A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#' is two folds.
#' @param weights An optional vector of length n containing sampling weights.
#' @param .bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .trim Determines the amount the density ratios should be trimmed.
#'  The default is 0.999, trimming the density ratios greater than the 0.999 percentile
#'  to the 0.999 percentile. A value of 1 indicates no trimming.
#' @param .SL_folds Integer. Controls the number of splits to be used for fitting
#'  the Super Learner. The default is 10.

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
#' \item{density_ratios}{An n x Tau matrix of the estimated density ratios.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' \item{outcome_type}{The outcome variable type.}
#'
#' @example inst/examples/tmle-ex.R
#' @export
lmtp_tmle <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                      cens = NULL, shift = NULL, shifted = NULL, k = Inf,
                      intervention_type = c("static", "dynamic", "mtp"),
                      outcome_type = c("binomial", "continuous", "survival"),
                      id = NULL, bounds = NULL, learners_outcome = "SL.glm",
                      learners_trt = "SL.glm", folds = 10, weights = NULL,
                      .bound = 1e-5, .trim = 0.999, .SL_folds = 10) {
  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = weights,
    bounds = bounds,
    bound = .bound
  )

  pb <- progressr::progressor(meta$tau*folds*2)

  ratios <- cf_r(
    meta$data, meta$shifted_data,
    meta$folds, meta$trt,
    cens, meta$risk, meta$tau,
    meta$node_list$trt, learners_trt,
    pb, meta$weights_r,
    match.arg(intervention_type),
    .SL_folds,
    .trim
  )

  estims <- cf_tmle(
    meta$data,
    meta$shifted_data,
    meta$folds,
    "xyz", # CHANGE THIS TO "*tmp_lmtp_outcome*",
    cens,
    meta$risk,
    meta$tau,
    meta$node_list$outcome,
    meta$outcome_type,
    meta$m, meta$m,
    ratio_tmle_ipw(ratios),
    learners_outcome,
    meta$weights,
    meta$weights_m,
    .SL_folds,
    pb
  )

  theta_dr(
    list(
      estimator = "TMLE",
      m = list(natural = estims$natural, shifted = estims$shifted),
      r = ratios$ratios,
      tau = meta$tau,
      folds = meta$folds,
      id = meta$id,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      weights = weights,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      weights_m = estims$sl_weights,
      weights_r = ratios$sl_weights,
      outcome_type = meta$outcome_type
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
#' @param data A data frame in wide format containing all necessary variables
#'  for the estimation problem.
#' @param trt A vector containing the column names of treatment variables ordered by time.
#' @param outcome The column name of the outcome variable. In the case of time-to-event
#'  analysis, a vector containing the columns names of intermediate outcome variables and the final
#'  outcome variable ordered by time. Only numeric values are allowed. If the outcome type
#'  is binary, data should be coded as 0 and 1.
#' @param baseline An optional vector containing the column names of baseline covariates to be
#'  included for adjustment at every time point.
#' @param time_vary A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as the number of time points of observation. If missingness in the outcome is
#'  present or if time-to-event outcome, must be provided.
#' @param shift A two argument function that specifies how treatment variables should be shifted.
#'  See examples for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param shifted An optional data frame, the same as in \code{data}, but modified according
#'  to the treatment policy of interest. If specified, \code{shift} is ignored.
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param intervention_type The intervetion type, should be one of \code{"static"},
#'   \code{"dynamic"}, \code{"mtp"}.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial, survival).
#' @param id An optional column name containing cluster level identifiers.
#' @param bounds An optional vector of the bounds for continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners_outcome A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param learners_trt A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#' is two folds.
#' @param weights An optional vector of length n containing sampling weights.
#' @param .bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .trim Determines the amount the density ratios should be trimmed.
#'  The default is 0.999, trimming the density ratios greater than the 0.999 percentile
#'  to the 0.999 percentile. A value of 1 indicates no trimming.
#' @param .SL_folds Integer. Controls the number of splits to be used for fitting
#'  the Super Learner. The default is 10.
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
#' \item{density_ratios}{An n x Tau matrix of the estimated density ratios.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' \item{outcome_type}{The outcome variable type.}
#' @export
#'
#' @example inst/examples/sdr-ex.R
lmtp_sdr <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                     cens = NULL, shift = NULL, shifted = NULL, k = Inf,
                     intervention_type = c("static", "dynamic", "mtp"),
                     outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL, learners_outcome = "SL.glm",
                     learners_trt = "SL.glm", folds = 10, weights = NULL,
                     .bound = 1e-5, .trim = 0.999, .SL_folds = 10) {
  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = NULL,
    bounds = bounds,
    bound = .bound
  )

  pb <- progressr::progressor(meta$tau*folds*2)

  ratios <- cf_r(
    meta$data, meta$shifted_data,
    meta$folds, meta$trt,
    cens, meta$risk, meta$tau,
    meta$node_list$trt, learners_trt,
    pb, meta$weights_r,
    match.arg(intervention_type),
    .SL_folds,
    .trim
  )

  estims <-
    cf_sdr(
      meta$data,
      meta$shifted_data,
      meta$folds,
      "xyz",
      cens,
      meta$risk,
      meta$tau,
      meta$node_list$outcome,
      meta$outcome_type,
      meta$m, meta$m,
      ratios$ratios,
      learners_outcome,
      meta$weights,
      meta$weights_m,
      .SL_folds,
      pb
    )

  theta_dr(
    list(
      estimator = "SDR",
      m = list(natural = estims$natural, shifted = estims$shifted),
      r = ratios$ratios,
      tau = meta$tau,
      folds = meta$folds,
      id = meta$id,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      weights = weights,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      weights_m = estims$sl_weights,
      weights_r = ratios$sl_weights,
      outcome_type = meta$outcome_type
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
#' @param data A data frame in wide format containing all necessary variables
#'  for the estimation problem.
#' @param trt A vector containing the column names of treatment variables ordered by time.
#' @param outcome The column name of the outcome variable. In the case of time-to-event
#'  analysis, a vector containing the columns names of intermediate outcome variables and the final
#'  outcome variable ordered by time. Only numeric values are allowed. If the outcome type
#'  is binary, data should be coded as 0 and 1.
#' @param baseline An optional vector containing the column names of baseline covariates to be
#'  included for adjustment at every time point.
#' @param time_vary A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as the number of time points of observation. If missingness in the outcome is
#'  present or if time-to-event outcome, must be provided.
#' @param shift A two argument function that specifies how treatment variables should be shifted.
#'  See examples for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param shifted An optional data frame, the same as in \code{data}, but modified according
#'  to the treatment policy of interest. If specified, \code{shift} is ignored.
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial, survival).
#' @param id An optional column name containing cluster level identifiers.
#' @param bounds An optional vector of the bounds for continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#'  is two folds.
#' @param weights An optional vector of length n containing sampling weights.
#' @param .bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .SL_folds Integer. Controls the number of splits to be used for fitting
#'  the Super Learner. The default is 10.
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
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{outcome_type}{The outcome variable type.}
#' @export
#'
#' @example inst/examples/sub-ex.R
lmtp_sub <- function(data, trt, outcome, baseline = NULL, time_vary = NULL, cens = NULL,
                     shift = NULL, shifted = NULL, k = Inf,
                     outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL, learners = "SL.glm", folds = 10,
                     weights = NULL, .bound = 1e-5, .SL_folds = 10) {
  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    learners_trt = NULL,
    learners_outcome = learners,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = NULL,
    bounds = bounds,
    bound = .bound
  )

  pb <- progressr::progressor(meta$tau*folds)

  estims <- cf_sub(
    meta$data,
    meta$shifted_data,
    meta$folds,
    "xyz",
    meta$node_list$outcome,
    cens,
    meta$risk,
    meta$tau,
    meta$outcome_type,
    learners,
    meta$m,
    pb,
    meta$weights_m,
    .SL_folds
  )

  theta_sub(
    eta = list(
      m = estims$m,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      folds = meta$folds,
      weights = weights,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      weights_m = estims$sl_weights,
      outcome_type = meta$outcome_type
    )
  )
}

#' LMTP IPW Estimator
#'
#' Inverse probability of treatment weighting estimator for the effects of traditional causal
#' effects and modified treatment policies for both point treatment and longitudinal data
#' with binary, continuous, or time-to-event outcomes. Supports binary, categorical, and continuous exposures.
#'
#' @param data A data frame in wide format containing all necessary variables
#'  for the estimation problem.
#' @param trt A vector containing the column names of treatment variables ordered by time.
#' @param outcome The column name of the outcome variable. In the case of time-to-event
#'  analysis, a vector containing the columns names of intermediate outcome variables and the final
#'  outcome variable ordered by time. Only numeric values are allowed. If the outcome type
#'  is binary, data should be coded as 0 and 1.
#' @param baseline An optional vector containing the column names of baseline covariates to be
#'  included for adjustment at every time point.
#' @param time_vary A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as the number of time points of observation. If missingness in the outcome is
#'  present or if time-to-event outcome, must be provided.
#' @param shift A two argument function that specifies how treatment variables should be shifted.
#'  See examples for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param shifted An optional data frame, the same as in \code{data}, but modified according
#'  to the treatment policy of interest. If specified, \code{shift} is ignored.
#' @param intervention_type The intervetion type, should be one of \code{"static"},
#'   \code{"dynamic"}, \code{"mtp"}.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial, survival).
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param id An optional column name containing cluster level identifiers.
#' @param learners A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#'  is two folds.
#' @param weights An optional vector of length n containing sampling weights.
#' @param .bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @param .trim Determines the amount the density ratios should be trimmed.
#'  The default is 0.999, trimming the density ratios greater than the 0.999 percentile
#'  to the 0.999 percentile. A value of 1 indicates no trimming.
#' @param .SL_folds Integer. Controls the number of splits to be used for fitting
#'  the Super Learner. The default is 10.
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
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' @export
#'
#' @example inst/examples/ipw-ex.R
lmtp_ipw <- function(data, trt, outcome, baseline = NULL, time_vary = NULL, cens = NULL,
                     shift = NULL, shifted = NULL,
                     intervention_type = c("static", "dynamic", "mtp"),
                     outcome_type = c("binomial", "continuous", "survival"),
                     k = Inf, id = NULL,
                     learners = "SL.glm", folds = 10, weights = NULL,
                     .bound = 1e-5, .trim = 0.999, .SL_folds = 10) {
  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    time_vary = time_vary,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    shifted = shifted,
    learners_trt = learners,
    learners_outcome = NULL,
    id = id,
    outcome_type = match.arg(outcome_type),
    V = folds,
    weights = NULL,
    bounds = NULL,
    bound = .bound
  )

  pb <- progressr::progressor(meta$tau*folds)

  ratios <- cf_r(
    meta$data, meta$shifted_data,
    meta$folds, meta$trt,
    cens, meta$risk, meta$tau,
    meta$node_list$trt, learners,
    pb, meta$weights_r,
    match.arg(intervention_type),
    .SL_folds,
    .trim
  )

  theta_ipw(
    eta = list(
      r = ratio_tmle_ipw(ratios),
      y = if (meta$survival) {
        convert_to_surv(data[[final_outcome(outcome)]])
      } else {
        data[[final_outcome(outcome)]]
      },
      folds = meta$folds,
      weights = weights,
      tau = meta$tau,
      shift = if (is.null(shifted)) deparse(substitute((shift))) else NULL,
      weights_r = ratios$sl_weights
    )
  )
}
