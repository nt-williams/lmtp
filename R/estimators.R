
#' LMTP Targeted Maximum Likelihood Estimator
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
#' @param pb Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # Estimating the effect of a point treatment
#' set.seed(6246)
#' n <- 500
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2))
#' Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
#' df <- data.frame(W, A, Y)
#' nodes <- list(c("W1", "W2"))
#' lmtp_tmle(df, "A", "Y", nodes, k = NULL, function(x) x + 2, "continuous",
#'           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
#'           learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))
lmtp_tmle <- function(data, A, Y, nodes, k = Inf, shift,
                      outcome_type = c("binomial", "continuous"),
                      bounds = NULL, learner_stack_Q = NULL,
                      learner_stack_g = NULL, pb = TRUE) {


  # setup
  n <- nrow(data)
  d <- shift_data(data, A, shift)
  t <- length(nodes)
  m <- matrix(nrow = n, ncol = t)
  ot <- match.arg(outcome_type)
  node_list <- create_node_list(A, nodes, k)
  scaled <- scale_y_continuous(data[, Y], ot, bounds)
  d[, "y_scaled"] <- data[, "y_scaled"] <- scaled$scaled
  pb_r <- check_pb(pb, t, "Estimating propensity")
  pb_m <- check_pb(pb, t, "Estimating regression")

  # propensity estimation
  r <- estimate_r(data, A, shift, t, node_list, learner_stack_g, pb_r)
  r <- use_dens_ratio(r, t, n, NULL, "tmle")

  # tmle
  m <- estimate_tmle(data = data,
                     shifted = d,
                     Y = "y_scaled",
                     node_list = node_list,
                     tau = t,
                     max = t,
                     outcome_type = ot,
                     m_natural = m,
                     m_shifted = m,
                     m_natural_initial = m,
                     m_shifted_initial = m,
                     r = r,
                     learner_stack = learner_stack_Q,
                     pb = pb_m)

  # estimates
  eta <- list(m = m,
              r = r,
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
#' @param pb Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' #' # Estimating the effect of a point treatment
#' set.seed(6246)
#' n <- 500
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2))
#' Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
#' df <- data.frame(W, A, Y)
#' nodes <- list(c("W1", "W2"))
#' lmtp_sdr(df, "A", "Y", nodes, k = NULL, function(x) x + 2, "continuous",
#'          learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
#'          learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))
lmtp_sdr <- function(data, A, Y, nodes, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learner_stack_Q = NULL,
                     learner_stack_g = NULL, pb = TRUE) {

  # setup
  n <- nrow(data)
  d <- shift_data(data, A, shift)
  t <- length(nodes)
  m <- matrix(nrow = n, ncol = t)
  ot <- match.arg(outcome_type)
  node_list <- create_node_list(A, nodes, k)
  scaled <- scale_y_continuous(data[, Y], ot, bounds)
  d[, "y_scaled"] <- data[, "y_scaled"] <- scaled$scaled
  pb_r <- check_pb(pb, t, "Estimating propensity")
  pb_m <- check_pb(pb, t, "Estimating regression")

  # propensity estimation
  r <- estimate_r(data, A, shift, t, node_list, learner_stack_g, pb_r)
  z <- use_dens_ratio(r, t, n, NULL, "eif")

  # sdr
  sdr <- estimate_sdr(data = data,
                      shifted = d,
                      Y = "y_scaled",
                      node_list = node_list,
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
#' @param pb Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # Estimating the effect of a point treatment
#' set.seed(6246)
#' n <- 500
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2))
#' Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
#' df <- data.frame(W, A, Y)
#' nodes <- list(c("W1", "W2"))
#' lmtp_sub(df, "A", "Y", nodes, k = NULL, function(x) x + 2, "continuous",
#'          learner_stack = sl3::make_learner(sl3::Lrnr_glm))
lmtp_sub <- function(data, A, Y, nodes, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, learner_stack = NULL, pb = TRUE) {

  # setup
  n <- nrow(data)
  d <- shift_data(data, A, shift)
  t <- length(nodes)
  m <- create_m(n, t, data[, Y])
  ot <- match.arg(outcome_type)
  node_list <- create_node_list(A, nodes, k)
  outcome_type = ifelse(ot == "continuous", "quasibinomial", "binomial")
  scaled <- scale_y_continuous(data[, Y], ot, bounds)
  d[, "y_scaled"] <- data[, "y_scaled"] <- scaled$scaled
  pb <- check_pb(pb, t, "Estimating regression")

  # sequential regression
  m <- estimate_sub(data = data,
                    shifted = d,
                    Y = "y_scaled",
                    node_list = node_list,
                    tau = t,
                    outcome_type = outcome_type,
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
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param shift A function that specifies how tratment variables should be shifted.
#' @param outcome_type Outcome variable type (i.e., continuous, binomial).
#' @param learner_stack An \code{sl3} learner stack for estimation of the
#'  exposure mechanism.
#' @param pb Should a progress bar be displayed? Default is \code{TRUE}.
#'
#' @return TODO
#' @export
#'
#' @examples
#' # Estimating the effect of a point treatment
#' set.seed(6246)
#' n <- 500
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2))
#' Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
#' df <- data.frame(W, A, Y)
#' nodes <- list(c("W1", "W2"))
#' lmtp_ipw(df, "A", "Y", nodes, k = NULL, function(x) x + 2, "continuous",
#'          learner_stack = sl3::make_learner(sl3::Lrnr_glm))
lmtp_ipw <- function(data, A, Y, nodes, k = Inf, shift,
                     outcome_type = c("binomial", "continuous"),
                     learner_stack = NULL, pb = TRUE) {

  # setup
  t <- length(nodes)
  n <- nrow(data)
  y <- data[, Y]
  node_list <- create_node_list(A, nodes, k)
  pb <- check_pb(pb, t, "Estimating propensity")

  # propensity estimation
  r <- estimate_r(data, A, shift, t, node_list, learner_stack, pb)
  r <- use_dens_ratio(r, t, n, NULL, "ipw")

  # estimates
  eta <- list(r = r,
              y = y,
              tau = t)

  out <- compute_theta(eta, "ipw")

  # returns
  return(out)

}
