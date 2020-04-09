
context("Fidelity of contrast estimates")
library(lmtp)
a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
stack <- sl3::make_learner(sl3::Lrnr_glm)

trt <- lmtp_tmle(sim_cens, a, "Y", nodes, baseline = NULL,
                 cens, k = 0, shift = function(x) x + 0.5,
                 outcome_type = "binomial",
                 learner_stack_Q = stack,
                 learner_stack_g = stack,
                 progress_bar = FALSE)

cntrl <- lmtp_tmle(sim_cens, a, "Y", nodes, baseline = NULL,
                   cens, k = 0, shift = NULL,
                   outcome_type = "binomial",
                   learner_stack_Q = stack,
                   learner_stack_g = stack,
                   progress_bar = FALSE)

lmtp_contrast(trt, cntrl, "additive")
lmtp_contrast(trt, cntrl, "relative risk")
lmtp_contrast(trt, cntrl, "odds ratio")
