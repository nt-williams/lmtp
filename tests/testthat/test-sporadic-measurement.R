context("Fidelity of estimates with sporadic outcome measurement")

n <- 1e4

set.seed(543634)

a1 <- rbinom(n, 1, 0.5)
s1 <- rbinom(n, 1, plogis(-2 - a1 * 0.25))
y1 <- rbinom(n, 1, plogis(-2 + a1 * 0.5))
a2 <- rbinom(n, 1, 0.5)
s2 <- rbinom(n, 1, plogis(-2 - a2 * 0.25 + y1 * 0.5))
y2 <- rbinom(n, 1, plogis(-2 + a2 * 0.5))
a3 <- rbinom(n, 1, 0.5)
y3 <- rbinom(n, 1, plogis(-2 + a3 * 0.5))

foo <- data.frame(a1 = a1, y1 = y1, a2 = a2, y2 = y2, a3 = a3, y3 = y3)
foo$y1[s1 == 1] <- NA_real_
foo$y2[s2 == 1] <- NA_real_

truth <- 0.119

curve <- lmtp_curve(foo, paste0("a", 1:3), paste0("y", 1:3), folds = 1, shift = static_binary_off)

# tests
test_that("Estimator fidelity with sporadic outcome and non-survival", {
  expect_equal(truth, curve$estimates[[2]]@x, tolerance = 0.01)
})

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
test_that("Estimator fidelity with sporadic outcome and survival", {
  expect_equal(truth, curve$estimates[[1]]@x, tolerance = 0.01)
})
