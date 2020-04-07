
#' LMTP Targeted Maximum Likelihood Estimator
#'
#' @param data A data frame.
#' @param A A vector of column names of treatment variables.
#' @param Y The column name of the outcome variable.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param cens A vector of column names of censoring indicators the same length as \code{A}.
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
lmtp_tmle <- function(data, A, Y, nodes, cens = NULL, k = Inf, shift,
                      outcome_type = c("binomial", "continuous"),
                      bounds = NULL, learner_stack_Q = NULL,
                      learner_stack_g = NULL, progress_bar = TRUE) {


  # setup
  n           <- nrow(data)
  shifted     <- shift_data(data, A, shift)
  t           <- length(nodes)
  m           <- matrix(nrow = n, ncol = t)
  ot          <- match.arg(outcome_type)
  node_list   <- create_node_list(A, nodes, k)
  scaled      <- scale_y_continuous(data[, Y], ot, bounds)
  shifted$xyz <- data$xyz <- scaled$scaled
  pb_r        <- check_pb(progress_bar, t, "Estimating propensity")
  pb_m        <- check_pb(progress_bar, t, "Estimating regression")

  # censoring
  cens_ratio <- estimate_c(data, cens, Y, t, node_list, learner_stack_g)

  # propensity estimation
  r <- estimate_r(data, A, cens, cens_ratio, shift, t, node_list, learner_stack_g, pb_r)
  z <- use_dens_ratio(r, t, n, NULL, "tml")

  # tmle engine
  m <- estimate_tmle(data = data,
                     shifted = shifted,
                     Y = "xyz",
                     node_list = node_list,
                     C = cens,
                     tau = t,
                     max = t,
                     outcome_type = ot,
                     m_natural = cbind(m, data[, Y]),
                     m_shifted = cbind(m, data[, Y]),
                     m_natural_initial = m,
                     m_shifted_initial = m,
                     r = z,
                     learner_stack = learner_stack_Q,
                     pb = pb_m)

  # estimates
  eta <- list(m = m,
              r = z,
              tau = t,
              outcome_type = ot,
              bounds = scaled$bounds)

  out <- compute_theta(eta, "tml")

  # returns
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
lmtp_sdr <- function(data, A, Y, nodes, cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learner_stack_Q = NULL,
                     learner_stack_g = NULL, progress_bar = TRUE) {

  # setup
  n           <- nrow(data)
  shifted     <- shift_data(data, A, shift)
  t           <- length(nodes)
  m           <- matrix(nrow = n, ncol = t)
  ot          <- match.arg(outcome_type)
  node_list   <- create_node_list(A, nodes, k)
  scaled      <- scale_y_continuous(data[, Y], ot, bounds)
  shifted$xyz <- data$xyz <- scaled$scaled
  pb_r        <- check_pb(progress_bar, t, "Estimating propensity")
  pb_m        <- check_pb(progress_bar, t, "Estimating regression")

  # censoring
  cens_ratio <- estimate_c(data, cens, Y, t, node_list, learner_stack_g)

  # propensity estimation
  r <- estimate_r(data, A, cens, cens_ratio, shift, t, node_list, learner_stack_g, pb_r)
  z <- use_dens_ratio(r, t, n, NULL, "eif")

  # sdr engine
  sdr <- estimate_sdr(data = data,
                      shifted = shifted,
                      Y = "xyz",
                      node_list = node_list,
                      C = cens,
                      tau = t,
                      max = t,
                      outcome_type = ot,
                      learner_stack = learner_stack_Q,
                      m_shifted = m,
                      m_natural = m,
                      r = r,
                      pb = pb_m)

  # estimates
  eta <- list(m = sdr,
              r = z,
              tau = t,
              outcome_type = ot,
              bounds = scaled$bounds)

  out <- compute_theta(eta, "sdr")

  # returns
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
lmtp_sub <- function(data, A, Y, nodes, cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learner_stack = NULL, progress_bar = TRUE) {

  # setup
  n           <- nrow(data)
  shifted     <- shift_data(data, A, shift)
  t           <- length(nodes)
  m           <- create_m(n, t, data[, Y])
  ot          <- match.arg(outcome_type)
  node_list   <- create_node_list(A, nodes, k)
  scaled      <- scale_y_continuous(data[, Y], ot, bounds)
  shifted$xyz <- data$xyz <- scaled$scaled
  pb          <- check_pb(progress_bar, t, "Estimating regression")

  # substitution engine
  m <- estimate_sub(data = data,
                    shifted = shifted,
                    Y = "xyz",
                    node_list = node_list,
                    C = cens,
                    tau = t,
                    outcome_type = ot,
                    m = m,
                    pb = pb)

  # estimates
  eta <- list(m = m,
              outcome_type = ot,
              bounds = scaled$bounds)

  out <- compute_theta(eta, "sub")

  # returns
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
#' @param cens A vector of column names of censoring indicators the same length as \code{A}.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how tratment variables should be shifted.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param learner_stack An \code{sl3} learner stack for estimation of the
#'  exposure mechanism.
#' @param progress_bar Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TO DO
lmtp_ipw <- function(data, A, Y, nodes, cens = NULL, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     learner_stack = NULL, progress_bar = TRUE) {

  # setup
  t         <- length(nodes)
  n         <- nrow(data)
  y         <- data[, Y]
  node_list <- create_node_list(A, nodes, k)
  pb        <- check_pb(progress_bar, t, "Estimating propensity")

  # censoring
  cens_ratio <- estimate_c(data, cens, Y, t, node_list, learner_stack)

  # propensity estimation
  r <- estimate_r(data, A, cens, cens_ratio, shift, t, node_list, learner_stack, pb)
  z <- use_dens_ratio(r, t, n, NULL, "ipw")

  # estimates
  eta <- list(r = z,
              y = y,
              tau = t)

  out <- compute_theta(eta, "ipw")

  # returns
  return(out)

}
