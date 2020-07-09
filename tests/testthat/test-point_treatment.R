
context("Fidelity of estimators for a point treatment")

# data generation
set.seed(429153)
n_obs <- 200
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
bs <- c("X1", "X2")
df <- data.frame(W, A, Y)
truth <- 0.76451

rule <- function(data, x) {
  data[[x]] + 0.5
}

# estimators
sub <-
  lmtp_sub(df, "A", "Y", baseline = bs, shift = rule,
           outcome_type = "binomial", folds = 2)

ipw <-
  lmtp_ipw(df, "A", "Y", baseline = bs, shift = rule, folds = 2)

tmle <-
  lmtp_tmle(df, "A", "Y", baseline = bs, shift = rule,
            outcome_type = "binomial", folds = 2)

sdr <-
  lmtp_sdr(df, "A", "Y", baseline = bs, shift = rule,
           outcome_type = "binomial", folds = 2)

# tests
test_that("point treatment fidelity", {
  expect_equal(truth, sub$theta, tolerance = 0.1)
  expect_equal(truth, ipw$theta, tolerance = 0.1)
  expect_equal(truth, tmle$theta, tolerance = 0.1)
  expect_equal(truth, sdr$theta, tolerance = 0.1)
})
