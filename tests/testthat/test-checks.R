
context("Various checks")

test_that("can't detect sl3", {
  verify_output(test_path("test-detect_sl3.txt"), {
    check_for_sl3(test = TRUE)
  })
})
