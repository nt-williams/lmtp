context("Fidelity of ltmle with categorical treatment")

set.seed(55)
n <- 1e4
W <- rbinom(n, 1, 0.5)

# Truth: potential outcomes under each level of a 3-level treatment
Y.0 <- rbinom(n, 1, plogis(-1.0 + 0.5 * W))
Y.1 <- rbinom(n, 1, plogis(-0.2 + 0.5 * W))
Y.2 <- rbinom(n, 1, plogis( 0.6 + 0.5 * W))
truth0 <- mean(Y.0)
truth1 <- mean(Y.1)
truth2 <- mean(Y.2)

# Observed data: treatment confounded by W
p <- cbind(plogis(-0.5 + W), plogis(0.3 - 0.4 * W), plogis(0.2 + 0.3 * W))
p <- p / rowSums(p)
A <- apply(p, 1, function(pr) sample(c("0", "1", "2"), 1, prob = pr))
Y <- ifelse(A == "0", Y.0, ifelse(A == "1", Y.1, Y.2))
tmp <- data.frame(W, A, Y)

tmle_cat <- sw(ltmle(tmp, "A", "Y", baseline = "W", folds = 1))

test_that("ltmle categorical treatment fidelity", {
  expect_equal(truth0, as.vector(tmle_cat$estimates$`0`@x), tolerance = 0.01)
  expect_equal(truth1, as.vector(tmle_cat$estimates$`1`@x), tolerance = 0.01)
  expect_equal(truth2, as.vector(tmle_cat$estimates$`2`@x), tolerance = 0.01)
})
