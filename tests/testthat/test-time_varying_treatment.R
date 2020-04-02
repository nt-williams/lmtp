
context("Fidelity of estimators for time-varying treatment")

# data generation, t = 2
set.seed(429)
n_obs <- 1000
L1 <- rbinom(n_obs, 1, 0.5)
A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))
df <- data.frame(L1, A1, L2, A2, Y)
a <- c("A1", "A2")
nodes <- list(c("L1"),
              c("L2"))
truth <- 0.88

# estimators
sub_glm <-
  lmtp_sub(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial", method = "glm")

sub_sl <-
  lmtp_sub(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial", method = "sl")

ipw <-
  lmtp_ipw(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner_stack(sl3::Lrnr_glm))

tmle <-
  lmtp_tmle(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
            outcome_type = "binomial",
            learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
            learner_stack_g = sl3::make_learner_stack(sl3::Lrnr_glm))

sdr <-
  lmtp_sdr(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
           learner_stack_g = sl3::make_learner_stack(sl3::Lrnr_glm))

# tests
test_that("time varying treatment fidelity", {
  expect_equal(truth, sub_glm, tolerance = 0.1)
  expect_equal(truth, sub_sl, tolerance = 0.1)
  expect_equal(truth, ipw, tolerance = 0.1)
  expect_equal(truth, tmle, tolerance = 0.1)
  expect_equal(truth, sdr, tolerance = 0.1)
})
