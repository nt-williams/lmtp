context("Fidelity of estimates with censoring")

A <- c("A1", "A2")
L <- list(c("L1"), c("L2"))
C <- c("C1", "C2")
truth <- 0.8

tmle <- sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = C, k = 0, shift = NULL, folds = 1))
sdr <- sw(lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = C, k = 0, shift = NULL, folds = 1))

# tests
test_that("estimator fidelity with censoring", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
})
