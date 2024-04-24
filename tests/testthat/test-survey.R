context("Survey weights")

set.seed(429153)

n <- 1e4
W1 <- rbinom(n, size = 1, prob = 0.5)
W2 <- rbinom(n, size = 1, prob = 0.65)
A <- rbinom(n, size = 1, prob = plogis(-0.4 + 0.2 * W2 + 0.15 * W1))
Y.1 <- rbinom(n, size = 1, prob = plogis(-1 + 1 - 0.1 * W1 + 0.3 * W2))
Y.0 <- rbinom(n, size = 1, prob = plogis(-1 + 0 - 0.1 * W1 + 0.3 * W2))

Y <- Y.1 * A + Y.0 * (1 - A)
tmp <- data.frame(W1, W2, A, Y, Y.1, Y.0)
truth <- mean(tmp$Y.1)

prob_S <- plogis(W1 * 0.5 + rnorm(n, mean = 0, sd = 1))
S <- rbinom(n, 1, prob_S)
tmp <- tmp[S == 1, ]
wts <- 1 / prob_S[S == 1]

sub <- lmtp_sub(tmp, "A", "Y", baseline = c("W1", "W2"), shift = static_binary_on,
                weights = wts, folds = 2)

ipw <- lmtp_ipw(tmp, "A", "Y", baseline = c("W1", "W2"), shift = static_binary_on,
                weights = wts, folds = 2)

tmle <- lmtp_tmle(tmp, "A", "Y", baseline = c("W1", "W2"), shift = static_binary_on,
                  weights = wts, folds = 2)

sdr <- lmtp_sdr(tmp, "A", "Y", baseline = c("W1", "W2"), shift = static_binary_on,
                weights = wts, folds = 2)


# tests
test_that("survey weight fidelity", {
  expect_equal(truth, sub$theta, tolerance = 0.025)
  expect_equal(truth, ipw$theta, tolerance = 0.025)
  expect_equal(truth, tmle$theta, tolerance = 0.025)
  expect_equal(truth, sdr$theta, tolerance = 0.025)
})
