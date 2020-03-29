
#' Parametric LMTP Substitution Estimator
#'
#' @param data data frame.
#' @param A a vector of column names of treatment variables.
#' @param Y column name of outcome variable.
#' @param history a list of length tau with the column names for treatment history at each time point.
#'   The list should be ordered following the time ordering of the model.
#' @param shift the amount to shift treatment variables.
#' @param family outcome variable family (i.e., continuous, binomial).
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(6246)
#' n <- 500
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2))
#' Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
#' df <- data.frame(W, A, Y)
#' history <- list(c("W1", "W2"))
#' lmtp_sub(df, "A", "Y", history, 2, "continuous")
lmtp_sub <- function(data, A, Y, history, shift,
                     outcome_type = c("binomial", "continuous"),
                     bounds = NULL) {

  # setup
  n <- nrow(data)
  tau <- length(history)
  node_list <- create_node_list(A, history)
  d <- shift_data(data, A, shift)
  m <- create_m(n, tau, data[, Y])
  ot <- match.arg(outcome_type)
  family <- ifelse(ot == "continuous", "gaussian", "binomial")
  scaled <- scale_y_continuous(data[, Y], ot, bounds)
  data[, "y_scaled"] <- scaled$scaled

  # iterative expectation estimation
  m <- estimate_m_glm(data, d, "y_scaled", node_list, tau, family, m)

  # estimation of theta
  theta <- compute_theta(m, "sub", ot, scaled$bounds)

  # returns
  return(theta)

}

lmtp_ipw <- function(data, A, Y, node_list, shift,
                     outcome_type = c("binomial", "continuous")) {

}
