
context("Fidelity of estimates with censoring")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.8

sub <-
  sw(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
              cens, k = 0, shift = NULL,
              outcome_type = "binomial",
              folds = 2))

ipw <-
  sw(lmtp_ipw(sim_cens, a, "Y", NULL, nodes,
              cens, k = 0, shift = NULL, folds = 2))

tmle <-
    sw(lmtp_tmle(sim_cens, a, "Y", nodes, baseline = NULL,
              cens, k = 0, shift = NULL,
              outcome_type = "binomial", folds = 2))

sdr <-
  sw(lmtp_sdr(sim_cens, a, "Y", nodes, baseline = NULL,
           cens, k = 0, shift = NULL,
           outcome_type = "binomial", folds = 2))

# tests
test_that("estimator fidelity with censoring", {
  expect_equal(truth, sub$theta, tolerance = 0.15)
  expect_equal(truth, ipw$theta, tolerance = 0.15)
  expect_equal(truth, tmle$theta, tolerance = 0.15)
  expect_equal(truth, sdr$theta, tolerance = 0.15)
})
