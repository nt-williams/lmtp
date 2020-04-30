
#' LMTP Targeted Maximum Likelihood Estimator
#'
#' @param data A data frame.
#' @param trt A vector of column names for treatment variables.
#' @param outcome The column name of the outcome variable.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint. If \code{k = Inf}, should be \code{NULL}
#'  and these variables should be added to the first index of \code{nodes}.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{A}.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how tratment variables should be shifted.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param bounds An optional vector of the bounds for continuous outcomes. If NULL
#'   the bounds will be taken as the minimum and maximum of the observed data.
#'   Ignored if outcome type is binary.
#' @param learners_outcome An \code{sl3} learner stack for estimation of the outcome
#'  regression.
#' @param learners_trt An \code{sl3} learner stack for estimation of the exposure
#'  mechanism.
#' @param folds The number of folds to be used for cross-validation.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_tmle <- function(data, trt, outcome, nodes, baseline = NULL,
                      cens = NULL, k = Inf, shift,
                      outcome_type = c("binomial", "continuous"),
                      bounds = NULL, learners_outcome = NULL,
                      learners_trt = NULL, folds = 10) {

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
    bounds = bounds
  )

  pb <- progressr::progressor(meta$tau*folds*2)

  # propensity --------------------------------------------------------------

  cens_rat <- cf_cens(
    data, meta$data, folds, cens, outcome,
    meta$tau, meta$node_list, learners_trt
  )

  dens_ratio <- ratio_dr(
    cf_r(meta$data, shift, folds, trt, cens, cens_rat,
         meta$tau, meta$node_list, learners_trt, pb),
    folds
  )

  # tmle --------------------------------------------------------------------

  estims <-
    cf_tmle(meta$data, meta$shifted_data, folds, "xyz", meta$node_list,
            cens, meta$tau, meta$outcome_type, meta$m,
            meta$m, dens_ratio, learners_outcome, pb)

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "tml",
    eta = list(
      m = estims,
      r = recombine_dens_ratio(dens_ratio),
      tau = meta$tau,
      folds = meta$folds,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      shift = deparse(substitute((shift)))
    ))

  return(out)
}

#' LMTP Sequential Doubly Robust Estimator
#'
#' @param data A data frame.
#' @param trt A vector of column names of treatment variables.
#' @param outcome The column name of the outcome variable.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint. If \code{k = Inf}, should be \code{NULL}
#'  and these variables should be added to the first index of \code{nodes}.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{A}.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how tratment variables should be shifted.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param bounds An optional vector of the bounds for continuous outcomes. If NULL
#'   the bounds will be taken as the minimum and maximum of the observed data.
#'   Ignored if outcome type is binary.
#' @param learners_outcome An \code{sl3} learner stack for estimation of the outcome
#'  regression.
#' @param learners_trt An \code{sl3} learner stack for estimation of the exposure
#'  mechanism.
#' @param folds The number of folds to be used for cross-validation.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_sdr <- function(data, trt, outcome, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learners_outcome = NULL,
                     learners_trt = NULL, folds = 10) {

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
    bounds = bounds
  )

  pb <- progressr::progressor(meta$tau*folds*2)

  # propensity --------------------------------------------------------------

  cens_rat <- cf_cens(data, meta$data, folds, cens, outcome,
                      meta$tau, meta$node_list, learners_trt)

  raw_ratio <- cf_r(meta$data, shift, folds, trt, cens, cens_rat,
                    meta$tau, meta$node_list, learners_trt, pb)

  # sdr ---------------------------------------------------------------------

  estims <-
    cf_sdr(meta$data, meta$shifted_data, folds, "xyz", meta$node_list,
           cens, meta$tau, meta$outcome_type, meta$m, meta$m, raw_ratio,
           learners_outcome, pb)

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "sdr",
    eta = list(
      m = estims,
      r = recombine_dens_ratio(ratio_dr(raw_ratio, folds)),
      tau = meta$tau,
      folds = meta$folds,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      shift = deparse(substitute((shift)))
    ))

  return(out)
}

#' LMTP Substitution Estimator
#'
#' @param data A data frame.
#' @param trt A vector of column names of treatment variables.
#' @param outcome The column name of the outcome variable.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint. If \code{k = Inf}, should be \code{NULL}
#'  and these variables should be added to the first index of \code{nodes}.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{A}.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how tratment variables should be shifted.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param bounds An optional vector of the bounds for continuous outcomes. If NULL
#'  the bounds will be taken as the minimum and maximum of the observed data.
#' @param learners An \code{sl3} learner stack for estimation of the outcome
#'  regression.
#' @param folds The number of folds to be used for cross-validation.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_sub <- function(data, trt, outcome, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learners = NULL, folds = 10) {

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
    bounds = bounds
  )

  pb <- progressr::progressor(meta$tau*folds)

  # substitution ------------------------------------------------------------

  estims <- cf_sub(meta$data, meta$shifted_data, folds, "xyz", meta$node_list,
                   cens, meta$tau, meta$outcome_type, learners, meta$m, pb)

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "sub",
    eta = list(
      m = estims,
      outcome_type = meta$outcome_type,
      bounds = meta$bounds,
      shift = deparse(substitute((shift)))
    ))

  return(out)

}

#' LMTP IPW Estimator
#'
#' @param data A data frame.
#' @param trt A vector of column names of treatment variables.
#' @param outcome The column name of the outcome variable.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint. If \code{k = Inf}, should be \code{NULL}
#'  and these variables should be added to the first index of \code{nodes}.
#' @param cens An optional vector of column names of censoring indicators the same
#'  length as \code{A}.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how tratment variables should be shifted.
#' @param learners An \code{sl3} learner stack for estimation of the
#'  exposure mechanism.
#' @param folds The number of folds to be used for cross-validation.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_ipw <- function(data, trt, outcome, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift, learners = NULL,
                     folds = 10) {

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
    bounds = NULL
  )

  pb <- progressr::progressor(meta$tau*folds)

  # propensity --------------------------------------------------------------

  cens_rat <- cf_cens(data, meta$data, folds, cens, outcome,
                      meta$tau, meta$node_list, learners)

  dens_ratio <-
    ratio_ipw(
      recombine_ipw(
        cf_r(meta$data, shift, folds, trt, cens, cens_rat,
             meta$tau, meta$node_list, learners, pb
        )
      )
    )

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    estimator = "ipw",
    eta = list(
      r = dens_ratio,
      y = data[[outcome]],
      folds = meta$folds,
      tau = meta$tau,
      shift = deparse(substitute((shift)))
    ))

  return(out)

}
