
context("Fidelity of estimators for time-varying treatment")

# data generation, t = 2
set.seed(429)
n_obs <- 200
L1 <- rbinom(n_obs, 1, 0.5)
A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))
df <- data.frame(L1, A1, L2, A2, Y)
a <- c("A1", "A2")
nodes <- list(c("L1"),
              c("L2"))
truth <- 0.88

rule <- function(data, x) {
  data[[x]] + 0.5
}

# estimators
sub <-
  lmtp_sub(df, a, "Y", time_vary = nodes, k = 0, shift = rule,
           outcome_type = "binomial",
           folds = 2)

ipw <-
  lmtp_ipw(df, a, "Y", time_vary = nodes, k = 0, shift = rule,
           folds = 2)

tmle <-
  lmtp_tmle(df, a, "Y", time_vary = nodes, cens = NULL, k = 0, shift = rule,
            outcome_type = "binomial",
            folds = 2)

sdr <-
  lmtp_sdr(df, a, "Y", time_vary = nodes, cens = NULL, k = 0, shift = rule,
           outcome_type = "binomial",
           folds = 2)

# tests
test_that("time varying treatment fidelity, t = 2", {
  expect_equal(truth, sub$theta, tolerance = 0.1)
  expect_equal(truth, ipw$theta, tolerance = 0.1)
  expect_equal(truth, tmle$theta, tolerance = 0.1)
  expect_equal(truth, sdr$theta, tolerance = 0.1)
})
