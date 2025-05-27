flip_sdr <- function(data, trt, outcome, baseline = NULL,
                     time_vary = NULL,
                     cens = NULL, compete = NULL,
                     overlap = FALSE, # if FALSE, use smooth trimming, else use overlap weights
                     trimming_threshold, #epsilon
                     smoothing_constant, # k
                     k = Inf,
                     outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL,
                     learners_outcome = "SL.glm",
                     learners_trt = "SL.glm",
                     folds = 10, weights = NULL,
                     control = lmtp_control(),
) {

}
