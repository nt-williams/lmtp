# Time-to-event analysis with a binary time-invariant exposure. Interested in
# the effect of treatment being given to all observations on the cumulative
# incidence of the outcome.
A <- "trt"
Y <- paste0("Y.", 1:6)
C <- paste0("C.", 0:5)
W <- c("W1", "W2")

lmtp_survival(sim_point_surv, A, Y, W, cens = C, folds = 1,
              shift = static_binary_on, estimator = "lmtp_tmle")
