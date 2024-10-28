context("Node list creation")

a <- c("A_1", "A_2", "A_3", "A_4")
bs <- c("W")
nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

foo = LmtpWideVars$new(NULL, nodes, a, NULL, "y", NULL, NULL, FALSE, 4, k = 0)

test_that("create proper node lists, t > 1", {
  verify_output(test_path("test-node-list.txt"), {
    a <- c("A_1", "A_2", "A_3", "A_4")
    bs <- c("W")
    nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

    # k = Inf
    foo <- LmtpWideVars$new(NULL, nodes, a, NULL, "y", NULL, NULL, FALSE, 4)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # k = 0
    foo <- LmtpWideVars$new(NULL, nodes, a, NULL, "y", NULL, NULL, FALSE, 4, k = 0)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # k = 1
    foo <- LmtpWideVars$new(NULL, nodes, a, NULL, "y", NULL, NULL, FALSE, 4, k = 1)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # k = 2
    foo <- LmtpWideVars$new(NULL, nodes, a, NULL, "y", NULL, NULL, FALSE, 4, k = 2)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # Markov
    foo <- LmtpWideVars$new(bs, nodes, a, NULL, "y", NULL, NULL, FALSE, 4, k = 0)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })
    create_node_list(a, 4, nodes, bs, k = 0)

    # Non-Markov with baseline
    foo <- LmtpWideVars$new(bs, nodes, a, NULL, "y", NULL, NULL, FALSE, 4)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })

    # NULL time-varying
    foo <- LmtpWideVars$new(bs, NULL, a, NULL, "y", NULL, NULL, FALSE, 4)
    lapply(1:4, function(t) {
      list(trt = foo$history("A", t),
           outcome = foo$history("L", t + 1))
    })
  })
})
