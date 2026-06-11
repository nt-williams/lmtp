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
tmle2 <- ltmle(tmp, "A", "Y", baseline = c("W1", "W2"), folds = 1)

# tests
test_that("point treatment fidelity", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
  expect_equal(truth, as.vector(tmle2$estimates$`1`@x), tolerance = 0.01)
})

set.seed(789)

n <- 1e4
W <- rbinom(n, 1, 0.5)

# truth: E[Y2 | do(A = 1)]
Y1.1 <- rbinom(n, 1, plogis(-3 + 0.5 * W + 1.5))
Y2.1 <- ifelse(Y1.1 == 1, 1, rbinom(n, 1, plogis(-2 + 0.3 * W + 1.2)))
truth_surv <- mean(Y2.1)

A <- rbinom(n, 1, plogis(-0.5 + 0.5 * W))
Y1 <- rbinom(n, 1, plogis(-3 + 0.5 * W + 1.5 * A))
Y2 <- ifelse(Y1 == 1, 1, rbinom(n, 1, plogis(-2 + 0.3 * W + 1.2 * A)))

tmp_surv <- data.frame(W, A, Y1, Y2)

tmle_surv <- sw(lmtp_tmle(tmp_surv, "A", c("Y1", "Y2"), baseline = "W",
                           shift = static_binary_on, outcome_type = "survival", folds = 1))
sdr_surv <- sw(lmtp_sdr(tmp_surv, "A", c("Y1", "Y2"), baseline = "W",
                         shift = static_binary_on, outcome_type = "survival", folds = 1))
tmle2 <- sw(ltmle(tmp_surv, "A", c("Y1", "Y2"), baseline = "W", outcome_type = "survival", folds = 1))

test_that("point treatment survival fidelity", {
  expect_equal(truth_surv, tmle_surv$estimate@x, tolerance = 0.01)
  expect_equal(truth_surv, sdr_surv$estimate@x, tolerance = 0.01)
})
