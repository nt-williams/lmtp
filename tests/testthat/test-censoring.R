
context("Fidelity of estimates with censoring")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.88

sub <-
  lmtp_sub(sim_cens, a, "Y", nodes, cens, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm),
           progress_bar = FALSE)

ipw <-
  lmtp_ipw(sim_cens, a, "Y", nodes, cens, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm),
           progress_bar = FALSE)

tmle <-
  lmtp_tmle(sim_cens, a, "Y", nodes, cens, k = 0, shift = function(x) x + 0.5,
            outcome_type = "binomial",
            learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
            learner_stack_g = sl3::make_learner(sl3::Lrnr_glm),
            progress_bar = FALSE)

sdr <-
  lmtp_sdr(sim_cens, a, "Y", nodes, cens, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
           learner_stack_g = sl3::make_learner(sl3::Lrnr_glm),
           progress_bar = FALSE)

# tests
test_that("estimator fidelity with censoring", {
  expect_equal(truth, sub$theta, tolerance = 0.1)
  expect_equal(truth, ipw$theta, tolerance = 0.1)
  expect_equal(truth, tmle$theta, tolerance = 0.1)
  expect_equal(truth, sdr$theta, tolerance = 0.1)
})
