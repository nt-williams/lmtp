context("Test helper functions")

ps <- matrix(nrow = 4, ncol = 2)
set.seed(3523)
ps[, 1] <- round(runif(4), 2)
ps[, 2] <- 1 - ps[, 1]
colnames(ps) <- c("0", "1")

set.seed(23)
a <- rbinom(4, 1, 0.5)

this_propensity(a, ps)

test_that("extracts correct propensity score for the observed level", {
  expect_equal(
    c(0.26, 0.07, 0.51, 0.06), this_propensity(a, ps)
  )
})
