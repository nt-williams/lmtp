
test_that("create proper node lists", {
  verify_output(test_path("test-node-list.txt"), {
    a <- c("A_1", "A_2", "A_3", "A_4")
    nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

    # k = null
    create_node_list(a, nodes, k = NULL)

    # k = 0
    create_node_list(a, nodes, k = 0)

    # k = 1
    create_node_list(a, nodes, k = 1)

    # k = 2
    create_node_list(a, nodes, k = 2)
  })
})
