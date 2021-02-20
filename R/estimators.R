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
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
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
#' @param return_all_ratios Logical. If \code{TRUE}, the non-cumulative product density
#'  ratios will be returned. The default is \code{FALSE}.
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
#' \item{raw_ratios}{An n x Tau matrix of the estimated non-cumulative product density ratios.
#'  \code{NULL} if \code{return_all_ratios = FALSE}.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' \item{outcome_type}{The outcome variable type.}
#'
#' @example inst/examples/tmle-ex.R
#' @export
lmtp_tmle <- function(data, trt, outcome, baseline = NULL,
                      time_vary = NULL, cens = NULL, shift, k = Inf,
                      outcome_type = c("binomial", "continuous", "survival"),
                      id = NULL, bounds = NULL, learners_outcome = "SL.glm",
                      learners_trt = "SL.glm", folds = 10, weights = NULL,
                      return_all_ratios = FALSE,
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

  ratios <- cf_r(meta$data, shift, folds, meta$trt, cens, meta$risk, meta$tau,
                 meta$node_list$trt, learners_trt, pb, meta$weights_r, .SL_folds)
  cumprod_ratios <- ratio_dr(ratios, folds, .trim)

  estims <-
    cf_tmle(meta$data, meta$shifted_data, folds, "xyz", meta$node_list$outcome,
            cens, meta$risk, meta$tau, meta$outcome_type, meta$m, meta$m,
            cumprod_ratios, learners_outcome, pb, meta$weights, meta$weights_m, .SL_folds)

  out <- compute_theta(
    estimator = "tml",
    eta = list(
      estimator = "TMLE",
      m = estims,
      r = recombine_dens_ratio(cumprod_ratios),
      raw_ratios = if (return_all_ratios) recombine_raw_ratio(ratios),
      tau = meta$tau,
      folds = meta$folds,
      id = meta$id,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      weights = weights,
      shift = deparse(substitute((shift))),
      weights_m = pluck_weights("m", estims),
      weights_r = pluck_weights("r", cumprod_ratios),
      outcome_type = meta$outcome_type
    ))

  return(out)
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
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
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
#' @param return_all_ratios Logical. If \code{TRUE}, the non-cumulative product density
#'  ratios will be returned. The default is \code{FALSE}.
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
#' \item{raw_ratios}{An n x Tau matrix of the estimated non-cumulative product density ratios.
#'  \code{NULL} if \code{return_all_ratios = FALSE}.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' \item{outcome_type}{The outcome variable type.}
#' @export
#'
#' @example inst/examples/sdr-ex.R
lmtp_sdr <- function(data, trt, outcome, baseline = NULL,
                     time_vary = NULL, cens = NULL, shift, k = Inf,
                     outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL, learners_outcome = "SL.glm",
                     learners_trt = "SL.glm", folds = 10, weights = NULL,
                     return_all_ratios = FALSE,
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

  ratios <- cf_r(meta$data, shift, folds, meta$trt, cens, meta$risk, meta$tau,
                 meta$node_list$trt, learners_trt, pb, meta$weights_r, .SL_folds)

  estims <-
    cf_sdr(meta$data, meta$shifted_data, folds, "xyz", meta$node_list$outcome,
           cens, meta$risk, meta$tau, meta$outcome_type, meta$m, meta$m,
           ratios, learners_outcome, pb, meta$weights_m, .trim, .SL_folds)

  out <- compute_theta(
    estimator = "sdr",
    eta = list(
      estimator = "SDR",
      m = estims,
      r = recombine_dens_ratio(ratio_dr(ratios, folds, .trim)),
      raw_ratios = if (return_all_ratios) recombine_raw_ratio(ratios),
      tau = meta$tau,
      folds = meta$folds,
      weights = weights,
      id = meta$id,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      shift = deparse(substitute((shift))),
      weights_m = pluck_weights("m", estims),
      weights_r = pluck_weights("r", ratios),
      outcome_type = meta$outcome_type
    ))

  return(out)
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
lmtp_sub <- function(data, trt, outcome, baseline = NULL,
                     time_vary = NULL, cens = NULL, shift, k = Inf,
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

  estims <- cf_sub(meta$data, meta$shifted_data, folds, "xyz",
                   meta$node_list$outcome,
                   cens, meta$risk, meta$tau, meta$outcome_type,
                   learners, meta$m, pb, meta$weights_m, .SL_folds)

  out <- compute_theta(
    estimator = "sub",
    eta = list(
      m = estims$m,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      folds = meta$folds,
      weights = weights,
      shift = deparse(substitute((shift))),
      weights_m = pluck_weights("m", estims),
      outcome_type = meta$outcome_type
    ))

  return(out)
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
#' @param return_all_ratios Logical. If \code{TRUE}, the non-cumulative product density
#'  ratios will be returned. The default is \code{FALSE}.
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
#' \item{raw_ratios}{An n x Tau matrix of the estimated non-cumulative product density ratios.
#'  \code{NULL} if \code{return_all_ratios = FALSE}.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' @export
#'
#' @example inst/examples/ipw-ex.R
lmtp_ipw <- function(data, trt, outcome, baseline = NULL,
                     time_vary = NULL, cens = NULL, k = Inf,
                     id = NULL, shift, outcome_type = c("binomial", "continuous", "survival"),
                     learners = "SL.glm", folds = 10, weights = NULL,
                     return_all_ratios = FALSE,
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

  ratios <- cf_r(meta$data, shift, folds, meta$trt, cens, meta$risk,
                 meta$tau, meta$node_list$trt, learners, pb, meta$weights_r,
                 .SL_folds)
  cumprod_ratios <- ratio_ipw(recombine_ipw(ratios), .trim)

  out <- compute_theta(
    estimator = "ipw",
    eta = list(
      r = cumprod_ratios$r,
      raw_ratios = if (return_all_ratios) recombine_raw_ratio(ratios),
      y = if (meta$survival) {
        convert_to_surv(data[[final_outcome(outcome)]])
      } else {
        data[[final_outcome(outcome)]]
      },
      folds = meta$folds,
      weights = weights,
      tau = meta$tau,
      shift = deparse(substitute((shift))),
      weights_r = cumprod_ratios$sl_weights
    ))

  return(out)
}
