
# point treatment
truth <- 22.97

sim_sub_glm <- function(n) {
  W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
  A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2)) + 2
  Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
  df <- data.frame(W, A, Y)
  lmtp_sub(df, "A", "Y", list(c("W1", "W2")), shift = function(x) x + 2,
           outcome_type = "continuous", method = "glm")
}

sim_sub_sl <- function(n) {
  W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
  A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2)) + 2
  Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
  df <- data.frame(W, A, Y)
  lmtp_sub(df, "A", "Y", list(c("W1", "W2")), shift = function(x) x + 2,
           outcome_type = "continuous", method = "sl",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm_fast))
}

sim_ipw <- function(n) {
  W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
  A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2)) + 2
  Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
  df <- data.frame(W, A, Y)
  lmtp_ipw(df, "A", "Y", list(c("W1", "W2")), shift = function(x) x + 2,
           outcome_type = "continuous",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm_fast))
}

sim_tmle <- function(n) {
  W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
  A <- rpois(n, lambda = exp(3 + .3*log(W$W1) - .2*exp(W$W1)*W$W2)) + 2
  Y <- rnorm(n, 1 + .5*A - .2*A*W$W2 + 2*A*tan(W$W1^2) - 2*W$W1*W$W2 + A*W$W1*W$W2, 1)
  df <- data.frame(W, A, Y)
  lmtp_tmle(df, "A", "Y", list(c("W1", "W2")), k = NULL, shift = function(x) x + 2,
           outcome_type = "continuous",
           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm_fast),
           learner_stack_g = sl3::make_learner(sl3::Lrnr_glm_fast))
}

set.seed(6246)
sub_test_glm <- mean(replicate(10, sim_sub_glm(20000)))

set.seed(4353)
sub_test_sl <- mean(replicate(10, sim_sub_sl(20000)))

set.seed(46243)
ipw_test <- mean(replicate(10, sim_ipw(20000)))

set.seed(753)
tmle_test <- mean(replicate(10, sim_tmle(1000)))

test_that("point treatment fidelity", {
  expect_equal(truth, sub_test_glm, tolerance = 0.1)
  expect_equal(truth, sub_test_sl, tolerance = 0.1)
  expect_equal(truth, ipw_test, tolerance = 0.1)
  expect_equal(truth, tmle_test, tolerance = 0.1)
})
