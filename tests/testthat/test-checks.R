
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
