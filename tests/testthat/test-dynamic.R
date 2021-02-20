
dynamic_vec <- function(data, trt) {
  # the function should either be vectorized or iterate over the rows
  if (trt == "A1") { # if the first time point set to 1
    return(rep(1, nrow(data)))
  } else { # a vectorized version
    (data[["L"]] > 0)*1 + (data[["L"]] <= 0)*0 # else return 1 or 0 based on L
  }
}

time_vary_on <- function(data, trt) {
  if (trt == "A1") return(rep(1, nrow(data)))
  else return(rep(0, nrow(data)))
}

# test adapted from the ltmle package vignette for a
# dynamic intervention example with censoring
rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
n <- 10000
W <- rnorm(n)
A1 <- rexpit(W)
C1 <- rexpit(0.6 * W - 0.5 * A1)
uncensored <- C1 == 1
L <- A2 <- C2 <- Y <- rep(NA, n)
L[uncensored] <- (0.3 * W[uncensored] + 0.2 * A1[uncensored] + rnorm(sum(uncensored)))
A2[uncensored] <- rexpit(W[uncensored] + A1[uncensored] + L[uncensored])
C2[uncensored] <- 1
C2[!uncensored] <- 0
Y[uncensored] <- rexpit(W[uncensored] - 0.6 * A1[uncensored] + L[uncensored] - 0.8 * A2[uncensored])
sim <- data.frame(W, A1, C1, L, A2, C2, Y)
a <- c("A1", "A2")
baseline <- "W"
cens <- c("C1", "C2")
nodes <- list(c(NULL), c("L"))

# static on or off at all time points
# truth = 0.308
tml.stc <- sw(lmtp_tmle(sim, a, "Y", baseline, nodes,
                     cens, shift = static_binary_on,
                     outcome_type = "binomial", folds = 2,
                     .SL_folds = 3))

# truth = 0.528
sdr.stc <- sw(lmtp_sdr(sim, a, "Y", baseline, nodes,
                    cens, shift = static_binary_off,
                    outcome_type = "binomial", folds = 2,
                    .SL_folds = 3))

# time varying dynamic
# truth = 0.433
tml.tv <- sw(lmtp_tmle(sim, a, "Y", baseline, nodes, cens, shift = time_vary_on,
                    outcome_type = "binomial", folds = 2, .SL_folds = 3))

sdr.tv <- sw(lmtp_sdr(sim, a, "Y", baseline, nodes, cens, shift = time_vary_on,
                   outcome_type = "binomial", folds = 2, .SL_folds = 3))

# time varying and covariate dynamic
# truth = 0.345
tml.dyn <- sw(lmtp_tmle(sim, a, "Y", baseline, nodes, cens, shift = dynamic_vec,
                     outcome_type = "binomial", folds = 2, .SL_folds = 3))

sdr.dyn <- sw(lmtp_sdr(sim, a, "Y", baseline, nodes, cens, shift = dynamic_vec,
                    outcome_type = "binomial", folds = 2, .SL_folds = 3))

# tests
test_that("Dynamic intervention fidelity", {
  expect_equal(0.308, tml.stc$theta, tolerance = 0.1)
  expect_equal(0.528, sdr.stc$theta, tolerance = 0.1)
  expect_equal(0.433, tml.tv$theta, tolerance = 0.1)
  expect_equal(0.433, sdr.tv$theta, tolerance = 0.1)
  expect_equal(0.345, tml.dyn$theta, tolerance = 0.1)
  expect_equal(0.345, sdr.dyn$theta, tolerance = 0.1)
})
