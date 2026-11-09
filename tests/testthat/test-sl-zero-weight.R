test_that("run_ensemble falls back to the best learner when NNLS gives every learner weight 0", {
  set.seed(1)
  n <- 84; tau <- 28; p <- 20
  d <- data.frame(id = 1:n, W1 = rnorm(n), W2 = rbinom(n, 1, 0.4))
  for (t in 0:(tau - 1)) {
    for (j in 1:p) d[[sprintf("L%d_%02d", j, t)]] <- if (j %% 3 == 0) rbinom(n, 1, 0.05) else rnorm(n)
    d[[sprintf("A_%02d", t)]] <- rbinom(n, 1, plogis(-1.5 + d[[sprintf("L1_%02d", t)]]))
    d[[sprintf("Y_%02d", t)]] <- rbinom(n, 1, 0.01)
  }
  A <- sprintf("A_%02d", 0:(tau - 1))
  Y <- sprintf("Y_%02d", 0:(tau - 1))
  L <- lapply(0:(tau - 1), function(t) sprintf("L%d_%02d", 1:p, t))
  for (t in 2:tau) d[[Y[t]]][d[[Y[t - 1]]] == 1] <- 1
  for (t in 2:tau) d[d[[Y[t - 1]]] == 1, c(A[t], L[[t]])] <- NA

  fit <- suppressWarnings(
    lmtp_tmle(d, trt = A, outcome = Y, baseline = c("W1", "W2"), time_vary = L,
              shifted = d, mtp = TRUE, k = 1, outcome_type = "survival", id = "id", folds = 1,
              learners_trt = "SL.glm", learners_outcome = "SL.glm",
              control = lmtp_control(.learners_trt_folds = 2, .learners_outcome_folds = 2))
  )
  # Natural course: the estimate must equal the observed survival at the horizon
  expect_equal(tidy(fit)$estimate[1], 1 - mean(d[[Y[tau]]]), tolerance = 1e-3)   # predictions are bounded at 1e-5
})
