
context("Fidelity of estimators for time-varying treatment")

# data generation, t = 2
set.seed(429)
n_obs <- 1000
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

# estimators
sub <-
  lmtp_sub(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial", learner_stack = sl3::make_learner(sl3::Lrnr_glm))

ipw <-
  lmtp_ipw(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner_stack(sl3::Lrnr_glm))

tmle <-
  lmtp_tmle(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
            outcome_type = "binomial",
            learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
            learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))

sdr <-
  lmtp_sdr(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
           learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))

# tests
test_that("time varying treatment fidelity, t = 2", {
  expect_equal(truth, sub$theta, tolerance = 0.1)
  expect_equal(truth, ipw$theta, tolerance = 0.1)
  expect_equal(truth, tmle$theta, tolerance = 0.1)
  expect_equal(truth, sdr$theta, tolerance = 0.1)
})

# t = 4
# a <- c("A_1", "A_2", "A_3", "A_4")
# nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
# d <- function(A) {
#   delta <- 1
#   return((A - delta) * (A - delta >= 0) + A * (A - delta < 0))
# }
#
# lmtp_sub(sim_t4, a, "Y", nodes, k = 0, shift = d, learner_stack = sl3::make_learner(sl3::Lrnr_glm))
#
# tau <- 4
# Ddag <- DAG.empty() +
#   node("L", t = 1, distr = "rcat.b1",
#        probs = c(0.5, 0.25, 0.25)) +
#   node("A", t = 1, distr = "rbinom", size = 5,
#        prob = (L[1] > 1) * 0.5 + (L[1] > 2) * 0.1) +
#   node("L", t = 2:tau, distr = "rbern",
#        prob = plogis(- 0.3 * L[t-1] + 0.5 * A[t-1])) +
#   node("A", t = 2:tau, distr = "rbinom", size = 5,
#        prob = plogis(1 / (1 + L[t] * 2 + A[t-1]))) +
#   node("Y", t = (tau + 1), distr = "rbern",
#        prob = plogis(1 / (1 - 1.2 * A[tau] + 0.3 * L[tau])), EFU = TRUE)
#
#
# datagen <- function(n, tau, D = Ddag) {
#
#   D <- set.DAG(D)
#   data <- sim(D, n = as.integer(n))
#
#   names(data)[substr(names(data), 1, 1) == 'Y'] <- 'Y'
#
#   return(data)
# }
#
# foo <- function(n) {
#   df <- suppressMessages(datagen(n, 4))
#   lmtp_ipw(df, a, "Y", nodes, k = Inf, shift = d,
#            outcome_type = "binomial",
#            learner_stack = sl3::make_learner_stack(sl3::Lrnr_glm))
# }
#
# foo(10000)
#
# . <- replicate(50, foo(1000), simplify = F)
#
# mean(map_dbl(., "theta"))
# hist(map_dbl(., "theta"))
#
# x <- map(., "low")
# y <- map(., "high")
# mean(map2_lgl(x, y, ~ between(0.48, .x, .y)))

