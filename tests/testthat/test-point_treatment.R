
context("Fidelity of estimators for a point treatment")

# data generation
set.seed(429153)
n_obs <- 1000
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
nodes <- list(c("X1", "X2"))
df <- data.frame(W, A, Y)
truth <- 0.76451

# estimators
sub_glm <-
  lmtp_sub(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial", method = "glm")

sub_sl <-
  lmtp_sub(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial", method = "sl")

ipw <-
  lmtp_ipw(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner_stack(sl3::Lrnr_glm))

tmle <-
  lmtp_tmle(df, "A", "Y", nodes, shift = function(x) x + 0.5,
            outcome_type = "binomial",
            learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
            learner_stack_g = sl3::make_learner_stack(sl3::Lrnr_glm))$theta

sdr <-
  lmtp_sdr(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
           learner_stack_g = sl3::make_learner_stack(sl3::Lrnr_glm))$theta

# tests
test_that("point treatment fidelity", {
  expect_equal(truth, sub_glm, tolerance = 0.1)
  expect_equal(truth, sub_sl, tolerance = 0.1)
  expect_equal(truth, ipw, tolerance = 0.1)
  expect_equal(truth, tmle, tolerance = 0.1)
  expect_equal(truth, sdr, tolerance = 0.1)
})
