
context("Use of contrast function")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")

rule <- function(data, x) {
  data[[x]] + 0.5
}

set.seed(58)

fit1 <-
  lmtp_tmle(sim_cens[1:500, ], a, "Y", nodes, baseline = NULL,
            cens, k = 0, shift = rule,
            outcome_type = "binomial",
            folds = 2)

set.seed(679)

fit0 <-
  lmtp_tmle(sim_cens[1:500, ], a, "Y", baseline = NULL, nodes,
            cens, k = 0, shift = NULL,
            outcome_type = "binomial",
            folds = 2)

test_that("contrast output is correct", {
  verify_output(test_path("test-contrast.txt"), {

    # 1 object vs scalar ref
    lmtp_contrast(fit1, ref = 0.787)

    # 2 objects vs scalar ref
    lmtp_contrast(fit1, fit0, ref = 0.787)

    # 1 object vs 1 object ref, additive
    lmtp_contrast(fit1, ref = fit0)

    # 1 object vs 1 object ref, rr
    lmtp_contrast(fit1, ref = fit0, type = "rr")

    # 1 object vs 1 object ref, or
    lmtp_contrast(fit1, ref = fit0, type = "or")
  })
})
