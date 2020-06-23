
context("Fidelity of estimates with censoring")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.88

rule <- function(data, x) {
  data[[x]] + 0.5
}

sub <-
<<<<<<< HEAD
    lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
             cens, k = 0, shift = rule,
=======
    lmtp_sub(sim_cens, a, "Y", baseline = NULL, nodes,
             cens, k = 0, shift = function(x) x + 0.5,
>>>>>>> devel
             outcome_type = "binomial",
             learners = sl3::make_learner(sl3::Lrnr_glm),
             folds = 2)

ipw <-
<<<<<<< HEAD
  lmtp_ipw(sim_cens, a, "Y", nodes, baseline = NULL,
           cens, k = 0, shift = rule,
=======
  lmtp_ipw(sim_cens, a, "Y", baseline = NULL, nodes,
           cens, k = 0, shift = function(x) x + 0.5,
>>>>>>> devel
           learners = sl3::make_learner(sl3::Lrnr_glm),
           folds = 2)

tmle <-
<<<<<<< HEAD
    lmtp_tmle(sim_cens, a, "Y", nodes, baseline = NULL,
              cens, k = 0, shift = rule,
=======
    lmtp_tmle(sim_cens, a, "Y", baseline = NULL, nodes,
              cens, k = 0, shift = function(x) x + 0.5,
>>>>>>> devel
              outcome_type = "binomial",
              learners_outcome = sl3::make_learner(sl3::Lrnr_glm),
              learners_trt = sl3::make_learner(sl3::Lrnr_glm),
              folds = 2)

sdr <-
<<<<<<< HEAD
  lmtp_sdr(sim_cens, a, "Y", nodes, baseline = NULL,
           cens, k = 0, shift = rule,
=======
  lmtp_sdr(sim_cens, a, "Y", baseline = NULL, nodes,
           cens, k = 0, shift = function(x) x + 0.5,
>>>>>>> devel
           outcome_type = "binomial",
           learners_outcome = sl3::make_learner(sl3::Lrnr_glm),
           learners_trt = sl3::make_learner(sl3::Lrnr_glm),
           folds = 2)

# tests
test_that("estimator fidelity with censoring", {
  expect_equal(truth, sub$theta, tolerance = 0.15)
  expect_equal(truth, ipw$theta, tolerance = 0.15)
  expect_equal(truth, tmle$theta, tolerance = 0.15)
  expect_equal(truth, sdr$theta, tolerance = 0.15)
})
