
#' Parametric LMTP Substitution Estimator
#'
#' @param data data frame.
#' @param A a vector of column names of treatment variables.
#' @param Y column name of outcome variable.
#' @param history a list of length tau with the column names for treatment history at each time point.
#'   The list should be ordered following the time ordering of the model.
#' @param shift a function specified how tratment variables should be shifted.
#' @param family outcome variable family (i.e., continuous, binomial).
#' @param bounds an optional vector of the bounds for continuous outcomes. If NULL
#'   the bounds will be taken as the minimum and maximum of the observed data.
#' @param method should estimation be through glm or Super Learner (using sl3). Default is glm.
#' @return
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
#' history <- list(c("W1", "W2"))
#' lmtp_sub(df, "A", "Y", history, function(x) x + 2, "continuous", method = "glm")
lmtp_sub <- function(data, A, Y, history, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL, method = c("glm", "sl")) {

  # setup
  n <- nrow(data)
  d <- shift_data(data, A, shift)
  t <- length(history)
  m <- create_m(n, t, data[, Y])
  ot <- match.arg(outcome_type)
  node_list <- create_node_list(A, history)
  family <- ifelse(ot == "continuous", "gaussian", "binomial")
  outcome_type = ifelse(ot == "continuous", "quasibinomial", "binomial")
  scaled <- scale_y_continuous(data[, Y], ot, bounds)
  method <- match.arg(method)
  d[, "y_scaled"] <- data[, "y_scaled"] <- scaled$scaled

  # sequential regression
  m <- switch(method,
              "glm" = estimate_m_glm(data = data, shifted = d, Y = "y_scaled", node_list = node_list, tau = t, family = family, m = m),
              "sl" = estimate_m_sl(data = data, shifted = d, Y = "y_scaled", node_list = node_list, tau = t, outcome_type = outcome_type, learners = NULL, m = m))

  # estimation of theta
  # TODO: options for the specific estimator should be wrapped up in a list
  # see lmtp_ipw for example
  theta <- compute_theta(m, "sub", ot, scaled$bounds, method)

  # returns
  return(theta)

}

#' lMTP IPW Estimator
#'
#' @param data
#' @param A
#' @param Y
#' @param history
#' @param shift
#' @param outcome_type
#' @param learners
#'
#' @return
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
#' history <- list(c("W1", "W2"))
#' lmtp_ipw(df, "A", "Y", history, function(x) x + 2, "binomial")
lmtp_ipw <- function(data, A, Y, history, shift,
                     outcome_type = c("binomial", "continuous"),
                     learners = NULL) {

  # setup
  t <- length(history)
  node_list <- create_node_list(A, history)

  # propensity estimation
  r <- estimate_r_sl(data, A, shift, t, node_list, learners)
  eta <- list(r = r$rn,
              y = data[, Y],
              tau = t)

  # estimation of theta
  theta <- compute_theta(eta, "ipw", NULL, NULL, NULL)

  # returns
  return(theta)

}
