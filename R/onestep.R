
#' LMPTP Substitution Estimator
#'
#' @param data data frame.
#' @param A column names of treatment variables.
#' @param Y column name of outcome variable.
#' @param node_list a node list that describes the time ordered structure of variables.
#' @param shift the amount to shift treatment variables.
#' @param family outcome variable family (i.e., continuous, binomial).
#' @param method TODO: the model to be used for estimation.
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(429153)
#' n_obs <- 1000
#' W <- replicate(2, rbinom(n_obs, 1, 0.5))
#' A <- rnorm(n_obs, mean = 2 * W, sd = 1)
#' Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
#' tl <- list(c("X1", "X2", "A"))
#' df <- data.frame(W, A, Y)
#' lmtp_sub(df, "A", "Y", tl, 0.5, "binomial")
lmtp_sub <- function(data, A, Y, node_list, shift,
                     outcome_type = c("binomial", "continuous")) {

  # setup
  n <- nrow(data)
  tau <- length(tl)
  d <- shift_data(data, A, shift)
  m <- create_m(n, tau, data[, Y])
  ot <- match.arg(outcome_type)
  family <- ifelse(ot == "continous", "gaussian", "binomial")

  # iterative expectation estimation
  m <- estimate_m_glm(data, d, tau, node_list, Y, m, family)

  # estimation of theta
  if (ot == "continuous") {
    theta <- mean(m[, 1])
  } else {
    theta <- mean(rexpit(m[, 1]))
  }

  return(theta)

}

lmtp_ipw <- function(data, A, Y, node_list, shift, family,
                     method = NULL) {

}
