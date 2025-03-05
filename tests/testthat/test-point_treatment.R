context("Fidelity of estimators for a point treatment")

set.seed(543523)

n <- 1e4
W1 <- rbinom(n, size = 1, prob = 0.5)
W2 <- rbinom(n, size = 1, prob = 0.65)
A <- rbinom(n, size = 1, prob = plogis(-0.4 + 0.2 * W2 + 0.15 * W1))
Y.1 <-rbinom(n, size = 1, prob = plogis(-1 + 1 - 0.1 * W1 + 0.3 * W2))
Y.0 <- rbinom(n, size = 1, prob = plogis(-1 + 0 - 0.1 * W1 + 0.3 * W2))

Y <- Y.1 * A + Y.0 * (1 - A)
tmp <- data.frame(W1, W2, A, Y, Y.1, Y.0)
truth <- mean(tmp$Y.1)

tmle <- lmtp_tmle(tmp, "A", "Y", baseline = c("W1", "W2"), shift = static_binary_on, folds = 1)
sdr <- lmtp_sdr(tmp, "A", "Y", baseline = c("W1", "W2"), shift = static_binary_on, folds = 1)

# tests
test_that("point treatment fidelity", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
})
