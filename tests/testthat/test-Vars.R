test_that("create proper node lists, t > 1", {
  verify_output(test_path("test-Vars.txt"), {
    a <- c("A_1", "A_2", "A_3", "A_4")
    bs <- c("W")
    nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

    # k = Inf
    foo <- LmtpVars$new(bs, nodes, a, NULL, NULL, "y", "binomial", 4)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # k = 0
    foo <- LmtpVars$new(bs, nodes, a, NULL, NULL, "y", "binomial", 4, k = 0)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # k = 1
    foo <- LmtpVars$new(bs, nodes, a, NULL, NULL, "y", "binomial", 4, k = 1)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # k = 2
    foo <- LmtpVars$new(bs, nodes, a, NULL, NULL, "y", "binomial", 4, k = 2)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # NULL time-varying
    foo <- LmtpVars$new(bs, NULL, a, NULL, NULL, "y", "binomial", 4)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # point-treatment survival
    foo <- LmtpVars$new(bs, NULL, "A", NULL, NULL, paste0("Y_", 1:4), "survival", 4)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })
  })
})

test_that("rename variable correctly", {
  a <- c("A_1", "A_2")
  bs <- c("W")
  nodes <- list(c("L_1_1", "L_1_2"), c("L_2_1", "L_2_2"))
  foo <- LmtpVars$new(bs, nodes, a, NULL, NULL, "y", "binomial", 2)

  expect_equal(as.vector(foo$rename("W")), "W")
  expect_equal(as.vector(foo$rename("A_1")), "..i..A_1")
  expect_equal(as.vector(foo$rename(c("L_1_1", "L_1_2"))), c("..i..L_1", "..i..L_2"))
  expect_equal(as.vector(foo$rename(c("L_2_1", "L_2_2"))), c("..i..L_1", "..i..L_2"))
  expect_equal(as.vector(foo$rename("y")), "..i..Y_1")
})

test_that("pulls proper variables for a given time", {
  a <- c("A_1", "A_2")
  bs <- c("W")
  nodes <- list(c("L_1_1", "L_1_2"), c("L_2_1", "L_2_2"))
  foo <- LmtpVars$new(bs, nodes, a, NULL, NULL, "y", "binomial", 2)

  expect_equal(foo$time(1), c("W", "L_1_1", "L_1_2", "A_1", "y"))
  expect_equal(foo$time(2), c("W", "L_2_1", "L_2_2", "A_2", "y"))
})
