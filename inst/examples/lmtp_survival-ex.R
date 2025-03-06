\donttest{
# Example 1.1
# Time-to-event analysis with a binary time-invariant exposure. Interested in
# the effect of treatment being given to all observations on the cumulative
# incidence of the outcome.
A <- "trt"
Y <- paste0("Y.", 1:6)
C <- paste0("C.", 0:5)
W <- c("W1", "W2")

curve <- lmtp_survival(sim_point_surv, A, Y, W, cens = C, folds = 1,
                       shift = static_binary_on, estimator = "lmtp_tmle")

tidy(curve)

# Example 1.2
# Time-to-event analysis with a binary time-invariant exposure and a competing-risk.
lmtp_survival(
  data = sim_competing_risks,
  trt = "A",
  cens = paste0("C", 1:5),
  compete = paste0("D", 1:5),
  baseline = paste0("W", 1:5),
  outcome = paste0("Y", 1:5),
  shift = static_binary_on,
  folds = 1
)
}
