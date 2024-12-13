test_that("pivot wide to long correctly", {
  verify_output(test_path("test-pivot.txt"), {
    foo <- LmtpVars$new(NULL,
                        list(c("L_1"), c("L_2"), c("L_3"), c("L_4")),
                        c("A_1", "A_2", "A_3", "A_4"),
                        NULL,
                        NULL,
                        "Y",
                        "binomial",
                        4)

    tmp <- sim_t4
    tmp$..i..lmtp_id <- tmp$ID
    tmp$ID <- NULL

    head(pivot(tmp, foo))

    foo <- LmtpVars$new(NULL,
                        list(c("L1"), c("L2")),
                        c("A1", "A2"),
                        c("C1", "C2"),
                        NULL,
                        "Y", "binomial", 2)

    tmp <- sim_cens
    tmp$..i..lmtp_id <- 1:nrow(tmp)

    head(pivot(tmp, foo))
  })
})
