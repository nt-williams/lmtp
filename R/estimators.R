
#' LMTP Targeted Maximum Likelihood Estimator
#'
#' @param data A data frame.
#' @param A A vector of column names for treatment variables.
#' @param Y The column name of the outcome variable.
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
#' @param learner_stack_Q An \code{sl3} learner stack for estimation of the outcome
#'  regression.
#' @param learner_stack_g An \code{sl3} learner stack for estimation of the exposure
#'  mechanism.
#' @param progress_bar Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_tmle <- function(data, A, Y, nodes, baseline = NULL,
                      cens = NULL, k = Inf, shift,
                      outcome_type = c("binomial", "continuous"),
                      bounds = NULL, learner_stack_Q = NULL,
                      learner_stack_g = NULL, progress_bar = TRUE) {

  # setup -------------------------------------------------------------------

  meta <- prepare_mbased(
    data = data,
    A = A,
    Y = Y,
    nodes = nodes,
    baseline = baseline,
    k = k,
    shift = shift,
    outcome_type = match.arg(outcome_type),
    bounds = bounds
  )

  # propensity --------------------------------------------------------------

  z <- use_dens_ratio(
    ratio = estimate_r(
      data = data,
      A = A,
      cens = cens,
      C = estimate_c(data, cens, Y, meta$tau, meta$node_list, learner_stack_g),
      shift = shift,
      tau = meta$tau,
      node_list = meta$node_list,
      learner_stack = learner_stack_g,
      pb = check_pb(progress_bar, meta$t, "Estimating propensity")
    ),
    tau = meta$tau,
    n = meta$n,
    max_tau = NULL,
    what_estim = "tml"
  )

  # tmle --------------------------------------------------------------------

  m <- estimate_tmle(
    data = meta$data,
    shifted = meta$shifted_data,
    Y = "xyz",
    node_list = meta$node_list,
    C = cens,
    tau = meta$tau,
    max = meta$tau,
    outcome_type = meta$outcome_type,
    m_natural = meta$m,
    m_shifted = meta$m,
    r = z,
    learner_stack = learner_stack_Q,
    pb = check_pb(progress_bar, meta$t, "Estimating regression")
  )

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    eta = list(
      m = m,
      r = z,
      tau = meta$tau,
      outcome_type = meta$outcome_type,
      bounds = meta$scale_meta$bounds,
      shift = deparse(substitute((shift)))
    ),
    estimator = "tml"
  )

  return(out)
}

#' LMTP Sequential Doubly Robust Estimator
#'
#' @param data A data frame.
#' @param A A vector of column names of treatment variables.
#' @param Y The column name of the outcome variable.
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
#' @param learner_stack_Q An \code{sl3} learner stack for estimation of the outcome
#'  regression.
#' @param learner_stack_g An \code{sl3} learner stack for estimation of the exposure
#'  mechanism.
#' @param progress_bar Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_sdr <- function(data, A, Y, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learner_stack_Q = NULL,
                     learner_stack_g = NULL, progress_bar = TRUE) {

  # setup -------------------------------------------------------------------

  meta <- prepare_mbased(
    data = data,
    A = A,
    Y = Y,
    nodes = nodes,
    baseline = baseline,
    k = k,
    shift = shift,
    outcome_type = match.arg(outcome_type),
    bounds = bounds
  )

  # propensity --------------------------------------------------------------

  r <- estimate_r(
    data = data,
    A = A,
    cens = cens,
    C = estimate_c(data, cens, Y, meta$tau, meta$node_list, learner_stack_g),
    shift = shift,
    tau = meta$tau,
    node_list = meta$node_list,
    learner_stack = learner_stack_g,
    pb = check_pb(progress_bar, meta$t, "Estimating propensity")
  )

  z <- use_dens_ratio(
    ratio = r,
    tau = meta$tau,
    n = meta$n,
    max_tau = NULL,
    what_estim = "eif"
  )

  # sdr ---------------------------------------------------------------------

  sdr <- estimate_sdr(
    data = meta$data,
    shifted = meta$shifted_data,
    Y = "xyz",
    node_list = meta$node_list,
    C = cens,
    tau = meta$tau,
    max = meta$tau,
    outcome_type = meta$outcome_type,
    learner_stack = learner_stack_Q,
    m_shifted = meta$m,
    m_natural = meta$m,
    r = r,
    pb = check_pb(progress_bar, meta$t, "Estimating regression")
  )

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    eta = list(
      m = sdr,
      r = z,
      tau = meta$tau,
      outcome_type = meta$outcome_type,
      bounds = meta$scale_meta$bounds,
      shift = deparse(substitute((shift)))
    ),
    estimator = "sdr"
  )

  return(out)
}

#' LMTP Substitution Estimator
#'
#' @param data A data frame.
#' @param A A vector of column names of treatment variables.
#' @param Y The column name of the outcome variable.
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
#' @param learner_stack An \code{sl3} learner stack for estimation of the outcome
#'  regression.
#' @param progress_bar Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_sub <- function(data, A, Y, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learner_stack = NULL, progress_bar = TRUE) {

  # setup -------------------------------------------------------------------

  meta <- prepare_mbased(
    data = data,
    A = A,
    Y = Y,
    nodes = nodes,
    baseline = baseline,
    k = k,
    shift = shift,
    outcome_type = match.arg(outcome_type),
    bounds = bounds
  )

  # substitution ------------------------------------------------------------

  m <- estimate_sub(
    data = meta$data,
    shifted = meta$shifted_data,
    Y = "xyz",
    node_list = meta$node_list,
    C = cens,
    tau = meta$tau,
    outcome_type = meta$outcome_type,
    m = meta$m,
    pb = check_pb(progress_bar, meta$t, "Estimating regression")
  )

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    eta = list(
      m = m,
      outcome_type = meta$outcome_type,
      bounds = meta$scale_meta$bounds,
      shift = deparse(substitute((shift)))
    ),
    estimator = "sub"
  )

  return(out)

}

#' LMTP IPW Estimator
#'
#' @param data A data frame.
#' @param A A vector of column names of treatment variables.
#' @param Y The column name of the outcome variable.
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
#' @param learner_stack An \code{sl3} learner stack for estimation of the
#'  exposure mechanism.
#' @param progress_bar Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_ipw <- function(data, A, Y, nodes, baseline = NULL,
                     cens = NULL, k = Inf, shift,
                     learner_stack = NULL, progress_bar = TRUE) {

  # setup -------------------------------------------------------------------

  meta <- prepare_rbased(data, A, Y, nodes, baseline, k, shift)

  # propensity --------------------------------------------------------------

  z <- use_dens_ratio(
    ratio = estimate_r(
      data = data,
      A = A,
      cens = cens,
      C = estimate_c(data, cens, Y, meta$t, meta$node_list, learner_stack),
      shift = shift,
      tau = meta$t,
      node_list = meta$node_list,
      learner_stack = learner_stack,
      pb = check_pb(progress_bar, meta$t, "Estimating propensity")
    ),
    tau = meta$tau,
    n = meta$n,
    max_tau = NULL,
    what_estim = "ipw"
  )

  # return estimates --------------------------------------------------------

  out <- compute_theta(
    eta = list(
      r = z,
      y = data[[Y]],
      tau = meta$tau,
      shift = deparse(substitute((shift)))
    ),
    estimator = "ipw"
  )

  return(out)

}
