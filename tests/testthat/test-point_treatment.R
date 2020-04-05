
context("Fidelity of estimators for a point treatment")

# data generation
set.seed(429153)
n_obs <- 1000
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
nodes <- list(c("X1", "X2"))
df <- data.frame(W, A, Y)
truth <- 0.76451

# estimators
sub <-
  lmtp_sub(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm))

ipw <-
  lmtp_ipw(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm))

tmle <-
  lmtp_tmle(df, "A", "Y", nodes, shift = function(x) x + 0.5,
            outcome_type = "binomial",
            learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
            learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))

sdr <-
  lmtp_sdr(df, "A", "Y", nodes, shift = function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
           learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))

# tests
test_that("point treatment fidelity", {
  expect_equal(truth, sub$theta, tolerance = 0.1)
  expect_equal(truth, ipw$theta, tolerance = 0.1)
  expect_equal(truth, tmle$theta, tolerance = 0.1)
  expect_equal(truth, sdr$theta, tolerance = 0.1)
})

# data generation

foo <- function(n) {
  n_obs <- n
  W <- replicate(2, rbinom(n_obs, 1, 0.5))
  A <- rnorm(n_obs, mean = 2 * W, sd = 1)
  cens <- rbinom(n_obs, 1, plogis(A * 3 + W + rnorm(n_obs, sd = 0.05)))
  # cens <- rep(1, n_obs)
  Y <- rep(NA, n_obs)
  Y[cens == 1] <- rbinom(sum(cens), 1, plogis(A[cens == 1] + W[cens == 1] + rnorm(sum(cens), mean = 0, sd = 1)))
  nodes <- list(c("X1", "X2"))
  df <- data.frame(W, A, cens, Y)
  ce <- estimate_c(df, "cens", "Y", 1, create_node_list("A", nodes, 0),
         sl3::make_learner(sl3::Lrnr_glm), matrix(nrow = n_obs, ncol = 1))
  r <- estimate_r(df, "A", ce, function(x) x + 0.5, 1,
                  create_node_list("A", nodes, 0),
                  sl3::make_learner(sl3::Lrnr_glm), NULL)

  y <- data.frame(r = use_dens_ratio(r, 1, n_obs, NULL, "ipw"), Y)
  mean(y$r * y$Y, na.rm = T)
}

foo(1000)

. <- replicate(100, foo(500))

set.seed(429153)
n_obs <- 500
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
cens <- rbinom(n_obs, 1, plogis(A * 3))
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
Y[cens == 0] <- NA
nodes <- list(c("X1", "X2"))
df <- data.frame(W, A, cens, Y)
truth <- 0.76451

c_fit <- (glm(cens ~ A + X1 + X2, data = df, family = "binomial"))
c_weights <- mean(cens) / predict(c_fit, type = "response")

r <- ratio_ite(estimate_r(df, "A", as.matrix(c_weights), function(x) x + 0.5, 1,
           create_node_list("A", nodes, 0),
           sl3::make_learner(sl3::Lrnr_glm), NULL), 1, n_obs)

y <- data.frame(r, Y)

mean(y$r * y$Y, na.rm = T)

# this actually works below!!
foo <- function(n) {
  n_obs <- n
  W <- replicate(2, rbinom(n_obs, 1, 0.5))
  A <- rnorm(n_obs, mean = 2 * W, sd = 1)
  cens <- rbinom(n_obs, 1, plogis(A * 3))
  Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
  Y[cens == 0] <- NA
  nodes <- list(c("X1", "X2"))
  df <- data.frame(W, A, cens, Y)
  truth <- 0.76451

  c_fit <- (glm(cens ~ A + X1 + X2, data = df, family = "binomial"))
  c_weights <- mean(cens) / predict(c_fit, type = "response")

  r <- ratio_ite(estimate_r(df, "A", as.matrix(c_weights), function(x) x + 0.5, 1,
                            create_node_list("A", nodes, 0),
                            sl3::make_learner(sl3::Lrnr_glm), NULL), 1, n_obs)

  y <- data.frame(r, Y)

  mean(y$r * y$Y, na.rm = T)
}

. <- replicate(500, foo(1000))


lmtp_ipw(df, "A", "Y", nodes, "cens", shift = function(x) x + 0.5,
         outcome_type = "binomial",
         learner_stack = sl3::make_learner(sl3::Lrnr_glm), pb = F)
