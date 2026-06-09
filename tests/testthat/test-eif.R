context("EIF properties")

set.seed(429)
n <- 500
W <- rbinom(n, 1, 0.5)
A <- rbinom(n, 1, plogis(-0.5 + 0.5 * W))
Y <- rbinom(n, 1, plogis(-1 + A + 0.5 * W))
dat <- data.frame(W, A, Y)

fit <- sw(lmtp_tmle(dat, "A", "Y", baseline = "W", shift = static_binary_on, folds = 5))

test_that("TMLE uncentered EIF mean equals point estimate under cross-fitting", {
  expect_equal(
    weighted.mean(fit$estimate@eif, fit$estimate@weights),
    fit$estimate@x,
    tolerance = 1e-10
  )
})
