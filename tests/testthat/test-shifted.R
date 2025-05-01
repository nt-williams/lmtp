context("Fidelity of estimates with shifted data supplied")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.88

sc <- shift_data(sim_cens, a, cens, function(data, trt) data[[trt]] + 0.5)

tmle <-
  sw(lmtp_tmle(sim_cens, a, "Y", nodes, baseline = NULL,
               cens, k = 0, shifted = sc,
               outcome_type = "binomial", folds = 2, mtp = TRUE))

sdr <-
  sw(lmtp_sdr(sim_cens, a, "Y", nodes, baseline = NULL,
              cens, k = 0, shifted = sc,
              outcome_type = "binomial", folds = 2, mtp = TRUE))

# tests
test_that("estimator fidelity with shifted data supplied", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
})
