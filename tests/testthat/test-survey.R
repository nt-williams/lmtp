context("Survey weights")

set.seed(429153)
n_obs <- 5000
W <- rbinom(n_obs, 1, 0.5)
prob_S <- plogis(W * 0.5 + rnorm(n_obs, mean = 0, sd = 1))
S <- rbinom(n_obs, 1, prob_S)
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))

df <- data.frame(W, A, Y)[S == 1, ]
truth <- 0.76562

wts <- 1 / prob_S[S == 1]

rule <- function(data, x) {
  data[[x]] + 0.5
}

sub <-
  lmtp_sub(df, "A", "Y", baseline = "W", shift = rule,
           outcome_type = "binomial", folds = 2, weights = wts)

ipw <-
  lmtp_ipw(df, "A", "Y", baseline = "W", shift = rule, folds = 2, weights = wts)

tmle <-
  lmtp_tmle(df, "A", "Y", baseline = "W", shift = rule,
            outcome_type = "binomial", folds = 2, weights = wts)

sdr <-
  lmtp_sdr(df, "A", "Y", baseline = "W", shift = rule,
           outcome_type = "binomial", folds = 2, weights = wts)

# tests
test_that("survey weight fidelity", {
  expect_equal(truth, sub$theta, tolerance = 0.1)
  expect_equal(truth, ipw$theta, tolerance = 0.1)
  expect_equal(truth, tmle$theta, tolerance = 0.1)
  expect_equal(truth, sdr$theta, tolerance = 0.1)
})
