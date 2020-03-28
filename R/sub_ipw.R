
#' Parametric LMPTP Substitution Estimator
#'
#' @param data data frame.
#' @param A column names of treatment variables.
#' @param Y column name of outcome variable.
#' @param node_list a node list that describes the time ordered structure of variables.
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
#' nl <- list(c("W1", "W2", "A"))
#' lmtp_sub(df, "A", "Y", nl, 2, "continuous")
lmtp_sub <- function(data, A, Y, node_list, shift,
                     outcome_type = c("binomial", "continuous")) {

  # setup
  n <- nrow(data)
  tau <- length(tl)
  d <- shift_data(data, A, shift)
  m <- create_m(n, tau, data[, Y])
  ot <- match.arg(outcome_type)
  family <- ifelse(ot == "continuous", "gaussian", "binomial")

  # iterative expectation estimation
  m <- estimate_m_glm(data, d, tau, node_list, Y, m, family)

  # estimation of theta
  # TODO: wrap this up into a separate helper function
  if (ot == "continuous") {
    theta <- mean(m[, 1])
  } else {
    theta <- mean(rexpit(m[, 1]))
  }

  return(theta)

}

lmtp_ipw <- function(data, A, Y, node_list, shift,
                     outcome_type = c("binomial", "continuous")) {

}
