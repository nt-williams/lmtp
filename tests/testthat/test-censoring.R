context("Fidelity of estimates with censoring")

A <- c("A1", "A2")
L <- list(c("L1"), c("L2"))
C <- c("C1", "C2")
truth <- 0.8

sub <- sw(lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = C, k = 0, shift = NULL, folds = 1))
ipw <- sw(lmtp_ipw(sim_cens, A, "Y", time_vary = L, cens = C, k = 0, shift = NULL, folds = 1))
tmle <- sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = C, k = 0, shift = NULL, folds = 1))
sdr <- sw(lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = C, k = 0, shift = NULL, folds = 1))

# tests
test_that("estimator fidelity with censoring", {
  expect_equal(truth, sub$theta, tolerance = 0.15)
  expect_equal(truth, ipw$theta, tolerance = 0.15)
  expect_equal(truth, tmle$theta, tolerance = 0.15)
  expect_equal(truth, sdr$theta, tolerance = 0.15)
})
