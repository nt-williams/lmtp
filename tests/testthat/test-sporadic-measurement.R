context("Fidelity of estimates with sporadic outcome measurement")

# Add sporadic missingness to the first time point
prob_sporadic <- with(sim_point_surv, plogis(-2 + W1*0.5 + W2 * 0.5))
sim_point_surv$Y.1[runif(nrow(sim_point_surv)) < prob_sporadic & sim_point_surv$C.1 == 1] <- NA_real_

A <- "trt"
Y <- paste0("Y.", 1:6)
C <- paste0("C.", 0:5)
W <- c("W1", "W2")

truth <- 0.99

curve <- lmtp_curve(sim_point_surv, A, Y, W, cens = C, folds = 1, shift = static_binary_on, outcome_type = "survival")

# tests
test_that("estimator fidelity with censoring", {
  expect_equal(truth, curve$estimates[[1]]@x, tolerance = 0.01)
})
