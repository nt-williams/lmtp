context("Fidelity of estimators with competing risks")

n <- 1e4

set.seed(4634)

suppressWarnings({
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W3 <- rnorm(n)

  D1 <- rep(0, n)
  Y1 <- rep(0, n)

  A1 <- rbinom(n, 1, plogis(0.3 * W1 - 0.2 * W2 + 0.1 * W3))
  C1 <- 1 - rbinom(n, 1, plogis(-4 + 0.1*W1 - 0.1*W2 + A1))

  D2 <- rbinom(n, 1, plogis(-3 + 0.2 * W1 - 0.1 * W2 + 0.1 * A1))
  Y2 <- ifelse(D2 == 1, 0, rbinom(n, 1, plogis(-3 + 0.2 * W1 - 0.1 * W2 + 0.3 * A1)))

  A2 <- ifelse(D2 == 1 | Y2 == 1, NA, rbinom(n, 1, plogis(0.2 * W1 - 0.1 * W2 + 0.3 * W3 + 0.4 * A1)))
  C2 <- ifelse(D2 == 1 | Y2 == 1, NA, 1 - rbinom(n, 1, plogis(-4 + 0.1* W1 - 0.1*W2 + 0.1*A1 + 0.2*A2)))

  D3 <- ifelse(D2 == 1, 1, ifelse(Y2 == 1, NA, rbinom(n, 1, plogis(-3 + 0.1 * W1 + 0.2 * W2 + 0.1 * A2))))
  Y3 <- ifelse(Y2 == 1, 1, ifelse(D3 == 1, 0, rbinom(n, 1, plogis(-3 + 0.3 * W1 - 0.2 * W2 + 0.4 * W3 + 0.5 * A2))))

  df <- data.frame(W1, W2, W3, A1, C1, D2, Y2, A2, C2, D3, Y3)

  df[c("D2", "Y2", "A2", "D3", "Y3")] <-
    lapply(df[c("D2", "Y2", "A2", "D3", "Y3")], function(x) ifelse(df$C1 == 0, NA, x))
  df$C2 <- ifelse(df$C1 == 0, 0, df$C2)
  df[c("D3", "Y3")] <- lapply(df[c("D3", "Y3")], function(x) ifelse(df$C2 == 0, NA, x))
  df[c("Y2", "Y3")] <- lapply(df[c("Y2", "Y3")], function(x) ifelse(df$D2 == 1, 0, x))
  df$D3 <- ifelse(df$D2 == 1, 1, df$D2)
  df$Y3 <- ifelse(df$Y2 == 1, 1, df$Y3)
})

tmle <- sw(lmtp_tmle(
  data = df,
  trt = c("A1", "A2"),
  compete = c("D2", "D3"),
  baseline = c("W1", "W2", "W3"),
  cens = c("C1", "C2"),
  outcome = c("Y2", "Y3"),
  shift = static_binary_on,
  outcome_type = "survival",
  folds = 1
))

sdr <- sw(lmtp_sdr(
  data = df,
  trt = c("A1", "A2"),
  compete = c("D2", "D3"),
  baseline = c("W1", "W2", "W3"),
  cens = c("C1", "C2"),
  outcome = c("Y2", "Y3"),
  shift = static_binary_on,
  outcome_type = "survival",
  folds = 1
))

# curve <- lmtp_curve(
#   data = df,
#   trt = c("A1", "A2"),
#   compete = c("D2", "D3"),
#   baseline = c("W1", "W2", "W3"),
#   cens = c("C1", "C2"),
#   outcome = c("Y2", "Y3"),
#   shift = static_binary_on,
#   outcome_type = "survival",
#   folds = 1
# )

truth <- 0.8688
test_that("estimator fidelity with competing risks", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
})
