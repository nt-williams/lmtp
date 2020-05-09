
context("Use of contrast function")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.88

set.seed(58)

fit1 <-
  lmtp_tmle(sim_cens[1:200, ], a, "Y", nodes, baseline = NULL,
            cens, k = 0, shift = function(x) x + 0.5,
            outcome_type = "binomial",
            learners_outcome = sl3::make_learner(sl3::Lrnr_glm),
            learners_trt = sl3::make_learner(sl3::Lrnr_glm),
            folds = 2)

set.seed(679)

fit0 <-
  lmtp_tmle(sim_cens[1:200, ], a, "Y", nodes, baseline = NULL,
            cens, k = 0, shift = NULL,
            outcome_type = "binomial",
            learners_outcome = sl3::make_learner(sl3::Lrnr_glm),
            learners_trt = sl3::make_learner(sl3::Lrnr_glm),
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
