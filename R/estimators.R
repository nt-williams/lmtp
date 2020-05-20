
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
#' @param nodes A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model. Must be provided, even
#'  if no time-varying covariates exist. In the case of a point-treatment, should be set
#'  to \code{list(c(NULL))}. If time-to-event with a time-invariant exposure,
#'  \code{nodes} should be the same length as the number of intermediate outcome variables
#'  with each index of the list similiarly set to \code{NULL}. See examples for demonstration.
#' @param baseline An optional vector of columns names of baseline covariates to be
#'  included for adjustment at every timepoint.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{nodes}. If missingness in the outcome is present, must be provided.
#' @param k An integer specifying how many previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how treatment variables should be shifted. See examples
#' for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param bounds An optional vector of the bounds for continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners_outcome An \code{sl3} learner stack for estimation of the outcome
#'  regression. If not specified, will default to an ensemble of an intercept only model
#'  and a GLM.
#' @param learners_trt An \code{sl3} learner stack for estimation of the exposure
#'  mechanism. If not specified, will default to an ensemble of an intercept only model
#'  and a GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#' is two folds.
#' @param bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @return A list of class \code{lmtp} containing the following components:
#'
#' \item{estimator}{The estimator used, in this case "TMLE".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{The estimated standard error of the LMTP effect.}
#' \item{low}{Lower bound of the 95% confidence interval of the LMTP effect.}
#' \item{high}{Upper bound of the 95% confidene interval of the LMTP effect.}
#' \item{eif}{The estimated, uncentered, influence function of the estimate.}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' \item{outcome_type}{The outcome variable type.}
#'
#' @example inst/examples/tmle-ex.R
#' @export
lmtp_tmle <- function(data, trt, outcome, nodes, baseline = NULL,
                      cens = NULL, k = Inf, shift,
                      outcome_type = c("binomial", "continuous"),
                      bounds = NULL, learners_outcome = NULL,
                      learners_trt = NULL, folds = 10, bound = 1e-5) {

  # setup -------------------------------------------------------------------

  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    nodes = nodes,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    outcome_type = match.arg(outcome_type),
    V = folds,
    bounds = bounds,
    bound = bound
  )

  pb <- progressr::progressor(meta$tau*folds*2)

  # propensity --------------------------------------------------------------

  dens_ratio <- ratio_dr(
    cf_r(meta$data, shift, folds, meta$trt, cens, meta$determ, meta$tau,
         meta$node_list, learners_trt, pb, meta$weights_r),
    folds
  )

  # tmle --------------------------------------------------------------------

  estims <-
    cf_tmle(meta$data, meta$shifted_data, folds, "xyz", meta$node_list,
            cens, meta$determ, meta$tau, meta$outcome_type, meta$m, meta$m,
            dens_ratio, learners_outcome, pb, meta$weights_m)

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "dr",
    eta = list(
      estimator = "TMLE",
      m = estims,
      r = recombine_dens_ratio(dens_ratio),
      tau = meta$tau,
      folds = meta$folds,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      shift = deparse(substitute((shift))),
      weights_m = pluck_weights("m", estims),
      weights_r = pluck_weights("r", dens_ratio),
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
#' @param nodes A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model. Must be provided, even
#'  if no time-varying covariates exist. In the case of a point-treatment, should be set
#'  to \code{list(c(NULL))}. If time-to-event with a time-invariant exposure,
#'  \code{nodes} should be the same length as the number of intermediate outcome variables
#'  with each index of the list similiarly set to \code{NULL}. See examples for demonstration.
#' @param baseline An optional vector of columns names of baseline covariates to be
#'  included for adjustment at every timepoint.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{nodes}. If missingness in the outcome is present, must be provided.
#' @param k An integer specifying how many previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how treatment variables should be shifted. See examples
#' for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param bounds An optional vector of the bounds for continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners_outcome An \code{sl3} learner stack for estimation of the outcome
#'  regression. If not specified, will default to an ensemble of an intercept only model
#'  and a GLM.
#' @param learners_trt An \code{sl3} learner stack for estimation of the exposure
#'  mechanism. If not specified, will default to an ensemble of an intercept only model
#'  and a GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#' is two folds.
#' @param bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#' @return A list of class \code{lmtp} containing the following components:
#'
#' \item{estimator}{The estimator used, in this case "SDR".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{The estimated standard error of the LMTP effect.}
#' \item{low}{Lower bound of the 95% confidence interval of the LMTP effect.}
#' \item{high}{Upper bound of the 95% confidene interval of the LMTP effect.}
#' \item{eif}{The estimated, uncentered, influence function of the estimate.}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' \item{outcome_type}{The outcome variable type.}
#' @export
#'
#' @example inst/examples/sdr-ex.R
lmtp_sdr <- function(data, trt, outcome, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learners_outcome = NULL,
                     learners_trt = NULL, folds = 10, bound = 1e-5) {

  # setup -------------------------------------------------------------------

  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    nodes = nodes,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    outcome_type = match.arg(outcome_type),
    V = folds,
    bounds = bounds,
    bound = bound
  )

  pb <- progressr::progressor(meta$tau*folds*2)

  # propensity --------------------------------------------------------------

  raw_ratio <- cf_r(meta$data, shift, folds, meta$trt, cens, meta$determ, meta$tau,
                    meta$node_list, learners_trt, pb, meta$weights_r)

  # sdr ---------------------------------------------------------------------

  estims <-
    cf_sdr(meta$data, meta$shifted_data, folds, "xyz", meta$node_list,
           cens, meta$determ, meta$tau, meta$outcome_type, meta$m, meta$m,
           raw_ratio, learners_outcome, pb, meta$weights_m)

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "dr",
    eta = list(
      estimator = "SDR",
      m = estims,
      r = recombine_dens_ratio(ratio_dr(raw_ratio, folds)),
      tau = meta$tau,
      folds = meta$folds,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      shift = deparse(substitute((shift))),
      weights_m = pluck_weights("m", estims),
      weights_r = pluck_weights("r", raw_ratio),
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
#' @param nodes A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model. Must be provided, even
#'  if no time-varying covariates exist. In the case of a point-treatment, should be set
#'  to \code{list(c(NULL))}. If time-to-event with a time-invariant exposure,
#'  \code{nodes} should be the same length as the number of intermediate outcome variables
#'  with each index of the list similiarly set to \code{NULL}. See examples for demonstration.
#' @param baseline An optional vector of columns names of baseline covariates to be
#'  included for adjustment at every timepoint.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{nodes}. If missingness in the outcome is present, must be provided.
#' @param k An integer specifying how many previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how treatment variables should be shifted. See examples
#'  for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param bounds An optional vector of the bounds for continuous outcomes. If \code{NULL},
#'  the bounds will be taken as the minimum and maximum of the observed data.
#'  Should be left as \code{NULL} if the outcome type is binary.
#' @param learners An \code{sl3} learner stack for estimation of the outcome
#'  regression. If not specified, will default to an ensemble of an intercept only model
#'  and a GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#'  is two folds.
#' @param bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#'
#' @return A list of class \code{lmtp} containing the following components:
#'
#' \item{estimator}{The estimator used, in this case "substitution".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{NA}
#' \item{low}{NA}
#' \item{high}{NA}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{weights_m}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the outcome regression.}
#' \item{outcome_type}{The outcome variable type.}
#' @export
#'
#' @example inst/examples/sub-ex.R
lmtp_sub <- function(data, trt, outcome, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learners = NULL, folds = 10, bound = 1e-5) {

  # setup -------------------------------------------------------------------

  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    nodes = nodes,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    outcome_type = match.arg(outcome_type),
    V = folds,
    bounds = bounds,
    bound = bound
  )

  pb <- progressr::progressor(meta$tau*folds)

  # substitution ------------------------------------------------------------

  estims <- cf_sub(meta$data, meta$shifted_data, folds, "xyz", meta$node_list,
                   cens, meta$determ, meta$tau, meta$outcome_type,
                   learners, meta$m, pb, meta$weights_m)

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "sub",
    eta = list(
      m = estims$m,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
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
#' @param nodes A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model. Must be provided, even
#'  if no time-varying covariates exist. In the case of a point-treatment, should be set
#'  to \code{list(c(NULL))}. If time-to-event with a time-invariant exposure,
#'  \code{nodes} should be the same length as the number of intermediate outcome variables
#'  with each index of the list similiarly set to \code{NULL}. See examples for demonstration.
#' @param baseline An optional vector of columns names of baseline covariates to be
#'  included for adjustment at every timepoint.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{nodes}. If missingness in the outcome is present, must be provided.
#' @param k An integer specifying how many previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how treatment variables should be shifted. See examples
#' for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param learners An \code{sl3} learner stack for estimation of the treatment mechanism.
#'  If not specified, will default to an ensemble of an intercept only model
#'  and a GLM.
#' @param folds The number of folds to be used for cross-fitting. Minimum allowable number
#'  is two folds.
#' @param bound Determines that maximum and minimum values (scaled) predictions
#'  will be bounded by. The default is 1e-5, bounding predictions by 1e-5 and 0.9999.
#'
#' @return A list of class \code{lmtp} containing the following components:
#'
#' \item{estimator}{The estimator used, in this case "IPW".}
#' \item{theta}{The estimated population LMTP effect.}
#' \item{standard_error}{NA}
#' \item{low}{NA}
#' \item{high}{NA}
#' \item{shift}{The shift function specifying the treatment policy of interest.}
#' \item{weights_r}{A list the same length as \code{folds}, containing the Super Learner
#'  ensemble weights at each time-point for each fold for the propensity.}
#' @export
#'
#' @example inst/examples/ipw-ex.R
lmtp_ipw <- function(data, trt, outcome, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift, learners = NULL,
                     folds = 10, bound = 1e-5) {

  # setup -------------------------------------------------------------------

  meta <- Meta$new(
    data = data,
    trt = trt,
    outcome = outcome,
    nodes = nodes,
    baseline = baseline,
    cens = cens,
    k = k,
    shift = shift,
    outcome_type = NULL,
    V = folds,
    bounds = NULL,
    bound = bound
  )

  pb <- progressr::progressor(meta$tau*folds)

  # propensity --------------------------------------------------------------

  dens_ratio <-
    ratio_ipw(
      recombine_ipw(
        cf_r(meta$data, shift, folds, meta$trt, cens, meta$deterministic,
             meta$tau, meta$node_list, learners, pb, meta$weights_r
        )
      )
    )

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "ipw",
    eta = list(
      r = dens_ratio$r,
      y = data[[final_outcome(outcome)]],
      folds = meta$folds,
      tau = meta$tau,
      shift = deparse(substitute((shift))),
      weights_r = dens_ratio$sl_weights
    ))

  return(out)

}
