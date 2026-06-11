context("Fidelity of estimators for time-varying treatment")

tmp <- sim_t4
a <- c("A_1", "A_2", "A_3", "A_4")
time_varying <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

for (i in a) {
  tmp[[i]] <- factor(tmp[[i]], levels = 0:5, ordered = FALSE)
}

d <- function(data, trt) {
  out <- list()
  a <- data[[trt]]
  for (i in 1:length(a)) {
    if (as.character(a[i]) %in% c("0", "1")) {
      out[[i]] <- as.character(a[i])
    } else {
      out[[i]] <- as.numeric(as.character(a[i])) - 1
    }
  }
  factor(unlist(out), levels = 0:5, ordered = FALSE)
}

truth <- 0.305

tmle <- sw(lmtp_tmle(tmp, a, "Y", time_vary = time_varying, shift = d, mtp = TRUE, folds = 1))
sdr <- sw(lmtp_sdr(tmp, a, "Y", time_vary = time_varying, shift = d, mtp = TRUE, folds = 1))

test_that("time varying treatment fidelity, t = 4", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
})

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
set.seed(234)
n <- 1e5
W <- rbinom(n, 1, 0.5)
A1 <- 0
L <- rexpit(0.3 * W + 0.2 * A1)
A2 <- 0
truth0 <- mean(rexpit(W - 0.6 * A1 + L - 0.8 * A2))

set.seed(234)
A1 <- 1
L <- rexpit(0.3 * W + 0.2 * A1)
A2 <- 1
truth1 <- mean(rexpit(W - 0.6 * A1 + L - 0.8 * A2))

set.seed(234)
n <- 1e4
W <- rbinom(n, 1, 0.5)
A1 <- rbinom(n, 1, 0.5)
L <- rexpit(0.3 * W + 0.2 * A1)
A2 <- rbinom(n, 1, 0.5)
Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
tmp <- data.frame(W, A1, L, A2, Y)

tmle2 <- ltmle(tmp, c("A1", "A2"), "Y", "W", list(NULL, "L"), folds = 1)

test_that("time varying treatment fidelity, ltmle", {
  expect_equal(truth0, as.vector(tmle2$estimates$`0`@x), tolerance = 0.01)
  expect_equal(truth1, as.vector(tmle2$estimates$`1`@x), tolerance = 0.01)
})
