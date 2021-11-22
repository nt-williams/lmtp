
context("Node list creation")

test_that("create proper node lists, t > 1", {
  verify_output(test_path("test-node-list.txt"), {

    a <- list(
      c("A_11", "A_21", "A_31", "A_41"),
      c("A_12", "A_22", "A_32", "A_42"),
      c("A_13", "A_23", "A_33", "A_43"),
      c("A_14", "A_24", "A_44", "A_44")
    )

    bs <- c("W")
    nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

    # k = NULL
    create_node_list(a, 4, nodes, baseline = NULL, k = NULL)

    # k = Inf
    create_node_list(a, 4, nodes, baseline = NULL, k = Inf)

    # k = 0
    create_node_list(a, 4, nodes, baseline = NULL, k = 0)

    # k = 1
    create_node_list(a, 4, nodes, baseline = NULL, k = 1)

    # k = 2
    create_node_list(a, 4, nodes, baseline = NULL, k = 2)

    # Markov
    create_node_list(a, 4, nodes, bs, k = 0)

    # Non-Markov with baseline
    create_node_list(a, 4, nodes, bs, k = Inf)

    # NULL time-varying
    create_node_list(a, 4, NULL, bs, k = Inf)
  })
})
