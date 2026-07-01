context("Fidelity of estimators for time-varying censoring with point-treatment")

# DGP: single point treatment A_1, censoring at 4 time points subject to
# time-varying confounding (L_2, L_3, L_4), single terminal outcome Y.
#
#   W -> A_1 -> C_1 -> L_2 -> C_2 -> L_3 -> C_3 -> L_4 -> C_4 -> Y
#              (A_1 also feeds every downstream L/C/Y directly)

generate <- function(n, intervene = NULL) {
  W   <- rnorm(n)
  A_1 <- if (!is.null(intervene)) rep(intervene, n) else rbinom(n, 1, plogis(0.3 * W))

  C_1 <- if (!is.null(intervene)) rep(1, n) else rbinom(n, 1, plogis(1.5 + 0.3 * W + 0.5 * A_1))

  L_2 <- rnorm(n, mean = 0.3 * W + 0.4 * A_1, sd = 1)
  C_2 <- if (!is.null(intervene)) rep(1, n) else rbinom(n, 1, plogis(1.5 + 0.3 * L_2 + 0.3 * A_1))

  L_3 <- rnorm(n, mean = 0.5 * L_2 + 0.3 * A_1, sd = 1)
  C_3 <- if (!is.null(intervene)) rep(1, n) else rbinom(n, 1, plogis(1.5 + 0.3 * L_3 + 0.3 * A_1))

  L_4 <- rnorm(n, mean = 0.5 * L_3 + 0.3 * A_1, sd = 1)
  C_4 <- if (!is.null(intervene)) rep(1, n) else rbinom(n, 1, plogis(1.5 + 0.3 * L_4 + 0.3 * A_1))

  Y <- rbinom(n, 1, plogis(-0.5 + 0.4 * L_4 + 0.8 * A_1 + 0.2 * W))

  df <- data.frame(W, A_1, C_1, L_2, C_2, L_3, C_3, L_4, C_4, Y)

  if (is.null(intervene)) {
    downstream_of_C1 <- c("L_2", "C_2", "L_3", "C_3", "L_4", "C_4", "Y")
    downstream_of_C2 <- c("L_3", "C_3", "L_4", "C_4", "Y")
    downstream_of_C3 <- c("L_4", "C_4", "Y")

    df[downstream_of_C1] <- lapply(df[downstream_of_C1], function(x) ifelse(df$C_1 == 0, NA, x))
    df[downstream_of_C2] <- lapply(df[downstream_of_C2], function(x) ifelse(df$C_2 == 0, NA, x))
    df[downstream_of_C3] <- lapply(df[downstream_of_C3], function(x) ifelse(df$C_3 == 0, NA, x))
    df$Y <- ifelse(df$C_4 == 0, NA, df$Y)
  }

  df
}

set.seed(1)
truth1 <- mean(generate(1e6, intervene = 1)$Y)
truth0 <- mean(generate(1e6, intervene = 0)$Y)

set.seed(2)
obs <- generate(1e4)

tmle_on <- sw(lmtp_tmle(
  obs, 
  trt = "A_1", 
  baseline = "W", 
  outcome = "Y",
  time_vary = list(NULL, "L_2", "L_3", "L_4"), 
  cens = c("C_1", "C_2", "C_3", "C_4"), 
  outcome_type = "binomial", 
  folds = 1, 
  shift = static_binary_on
))

sdr_off <- sw(lmtp_sdr(
  obs, 
  trt = "A_1", 
  baseline = "W", 
  outcome = "Y",
  time_vary = list(NULL, "L_2", "L_3", "L_4"), 
  cens = c("C_1", "C_2", "C_3", "C_4"), 
  outcome_type = "binomial", 
  folds = 1, 
  shift = static_binary_off
))

lt <- sw(ltmle(
  obs,
  trt = "A_1",
  outcome = "Y",
  baseline = "W",
  time_vary = list(NULL, "L_2", "L_3", "L_4"),
  cens = c("C_1", "C_2", "C_3", "C_4"),
  outcome_type = "binomial",
  folds = 1
))

test_that("time varying censoring fidelty", {
  expect_equal(truth1, tmle_on$estimate@x, tolerance = 0.01)
  expect_equal(truth0, sdr_off$estimate@x, tolerance = 0.01)
  expect_equal(truth1, as.vector(lt$estimates$`1`@x), tolerance = 0.01)
  expect_equal(truth0, as.vector(lt$estimates$`0`@x), tolerance = 0.01)
})
