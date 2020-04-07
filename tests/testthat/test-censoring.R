# data generation
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

# data generation, t = 2
set.seed(429)
n_obs <- 1000
L1 <- rbinom(n_obs, 1, 0.5)
A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
C1 <- rbinom(n_obs, 1, plogis(A1 * 3))
L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
C2 <- rbinom(n_obs, 1, plogis(A2 * 3))
Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))

L2[C1 == 0] <- NA
A2[C1 == 0] <- NA
C2[C1 == 0] <- 0
Y[C1 == 0] <- NA
Y[C2 == 0] <- NA

df <- data.frame(L1, A1, C1, L2, A2, C2, Y)
a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.88

# substitution
estimate_sub(df, shift_data(df, a, function(x) x + 0.5),
             "Y", create_node_list(a, nodes, 0), cens,
             2, "binomial", learner_stack = sl3::make_learner(sl3::Lrnr_glm),
             matrix(nrow = 688, ncol = 2), NULL)

foo <- function(n) {
  n_obs <- n
  L1 <- rbinom(n_obs, 1, 0.5)
  A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
  C1 <- rbinom(n_obs, 1, plogis(A1 * 3))
  L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
  A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
  C2 <- rbinom(n_obs, 1, plogis(A2 * 3))
  Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))

  L2[C1 == 0] <- NA
  A2[C1 == 0] <- NA
  C2[C1 == 0] <- 0
  Y[C1 == 0] <- NA
  Y[C2 == 0] <- NA

  df <- data.frame(L1, A1, C1, L2, A2, C2, Y)
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  # substitution
  mean(estimate_sub(df, shift_data(df, a, function(x) x + 0.5),
               "Y", create_node_list(a, nodes, 0), cens,
               2, "binomial", learner_stack = sl3::make_learner(sl3::Lrnr_glm),
               matrix(nrow = sum(complete.cases(df)), ncol = 2), NULL)[, 1])
}

. <- replicate(100, foo(1000))


# this actually works below for IPW
foo <- function(n) {
  n_obs <- n
  W <- replicate(2, rbinom(n_obs, 1, 0.5))
  A <- rnorm(n_obs, mean = 2 * W, sd = 1)
  cens <- rbinom(n_obs, 1, plogis(A * 3))
  Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
  Y[cens == 0] <- NA
  nodes <- list(c("X1", "X2"))
  df <- data.frame(W, A, cens, Y)
  lmtp_ipw(df, "A", "Y", nodes, "cens", k = 0, function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm),
           pb = FALSE)$theta
}

foo <- function(n) {
  set.seed(429)
  n_obs <- n
  L1 <- rbinom(n_obs, 1, 0.5)
  A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
  C1 <- rbinom(n_obs, 1, 1 - 0.1*L1)
  L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
  A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
  C2 <- rbinom(n_obs, 1, 1 - 0.1*L2)
  Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))

  L2[C1 == 0] <- NA
  A2[C1 == 0] <- NA
  C2[C1 == 0] <- 0
  Y[C1 == 0] <- NA
  Y[C2 == 0] <- NA

  df <- data.frame(L1, A1, C1, L2, A2, C2, Y)
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  lmtp_ipw(df, a, "Y", nodes, cens, k = 0, function(x) x + 0.5,
           outcome_type = "binomial",
           learner_stack = sl3::make_learner(sl3::Lrnr_glm),
           pb = FALSE)$theta
}

. <- replicate(50, foo(500))

# tmle
n_obs <- 1000
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
cens <- rbinom(n_obs, 1, plogis(A * 3))
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
Y[cens == 0] <- NA
nodes <- list(c("X1", "X2"))
df <- data.frame(W, A, cens, Y)

lmtp_tmle(df, "A", "Y", nodes, "cens", 0, function(x) x + 0.5,
          "binomial", NULL, sl3::make_learner(sl3::Lrnr_glm),
          sl3::make_learner(sl3::Lrnr_glm))$theta

foo <- function(n) {
  n_obs <- n
  W <- replicate(2, rbinom(n_obs, 1, 0.5))
  A <- rnorm(n_obs, mean = 2 * W, sd = 1)
  cens <- rbinom(n_obs, 1, plogis(A * 3))
  Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
  Y[cens == 0] <- NA
  nodes <- list(c("X1", "X2"))
  df <- data.frame(W, A, cens, Y)

  lmtp_tmle(df, "A", "Y", nodes, "cens", 0, function(x) x + 0.5,
            "binomial", NULL, sl3::make_learner(sl3::Lrnr_glm),
            sl3::make_learner(sl3::Lrnr_glm))$theta
}

. <- replicate(50, foo(1000))

foo <- function(n) {
  n_obs <- n
  L1 <- rbinom(n_obs, 1, 0.5)
  A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
  C1 <- rbinom(n_obs, 1, 1 - 0.1*L1)
  L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
  A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
  C2 <- rbinom(n_obs, 1, 1 - 0.1*L2)
  Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))

  L2[C1 == 0] <- NA
  A2[C1 == 0] <- NA
  C2[C1 == 0] <- 0
  Y[C1 == 0] <- NA
  Y[C2 == 0] <- NA

  df <- data.frame(L1, A1, C1, L2, A2, C2, Y)
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  lmtp_tmle(df, a, "Y", nodes, cens, 0, function(x) x + 0.5,
            "binomial", NULL, sl3::make_learner(sl3::Lrnr_glm),
            sl3::make_learner(sl3::Lrnr_glm))
}

. <- replicate(200, foo(1000), simplify = F)


# sdr
foo <- function(n) {
  n_obs <- n
  L1 <- rbinom(n_obs, 1, 0.5)
  A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
  C1 <- rbinom(n_obs, 1, 1 - 0.1*L1)
  L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
  A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
  C2 <- rbinom(n_obs, 1, 1 - 0.1*L2)
  Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))

  L2[C1 == 0] <- NA
  A2[C1 == 0] <- NA
  C2[C1 == 0] <- 0
  Y[C1 == 0] <- NA
  Y[C2 == 0] <- NA

  df <- data.frame(L1, A1, C1, L2, A2, C2, Y)
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  lmtp_sdr(df, a, "Y", nodes, cens, 0, function(x) x + 0.5,
           "binomial", NULL, sl3::make_learner(sl3::Lrnr_glm),
           sl3::make_learner(sl3::Lrnr_glm), progress_bar = F)
}

. <- replicate(100, foo(1000), simplify = F)
