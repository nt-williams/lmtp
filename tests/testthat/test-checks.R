
context("Various checks")

test_that("detects xyz", {
  df <- data.frame(xyz = 1:5)
  expect_error(check_scaled_conflict(df))
})

test_that("detects incorrect folds", {
  expect_error(check_folds(1))
})

test_that("variables dont exist", {
  a <- c("A", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  expect_error(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
                        cens, k = 1, shift = function(x) x + 0.5))
})

test_that("time_vary is a list", {
  a <- c("A1", "A2")
  nodes <- c("L1", "L2")
  cens <- c("C1", "C2")
  expect_error(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
                        cens, k = 1, shift = function(x) x + 0.5))
})

test_that("variable length mismatch", {
  a <- c("A1")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  expect_error(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
                        cens, k = 1, shift = function(x) x + 0.5))

  a <- c("A1", "A2")
  nodes <- list(c("L1"))
  expect_error(lmtp_sub(sim_cens[complete.cases(sim_cens), ], a, "Y",
                        nodes, baseline = NULL, k = 1, shift = function(x) x + 0.5))

})
