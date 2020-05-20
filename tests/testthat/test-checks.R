
context("Various checks")

test_that("can't detect sl3", {
  verify_output(test_path("test-detect_sl3.txt"), {
    check_for_sl3(test = TRUE)
  })
})

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
  truth <- 0.88
  expect_error(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
                        cens, k = 1, shift = function(x) x + 0.5))
})
