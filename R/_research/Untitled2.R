data("iptwExWide", package = "twang")

trt <- paste0("tx", 1:3)
baseline <- c("gender", "age")
time_vary <- list(c("use0"), c("use1"), c("use2"))
cens <- NULL
compete <- NULL
id <- NULL
outcome <- "outcome"

variable_names <- c(unlist(trt), outcome, unlist(time_vary), baseline, cens, compete, id)

task <- LmtpTask$new(
  data = iptwExWide,
  shifted = make_shifted(iptwExWide[, variable_names], trt, cens, static_binary_on, NULL),
  A = trt,
  Y = outcome,
  L = time_vary,
  W = baseline,
  C = cens,
  D = compete,
  k = Inf, id = id,
  outcome_type = "continuous",
  folds = 10,
  weights = NULL,
  bounds = NULL
)

progressr::with_progress({
  progress_bar <- progressr::progressor(task$time_horizon*10*2)
  propensity <- cf_propensity_score(task, "SL.glm", lmtp_control(), progress_bar)
})

estimate_tmle_delay(task, 1, propensity$propensity_score, "SL.glm", lmtp_control(), NULL)


delay_weights(task$natural$tx3, propensity$propensity_score, "1", 3, 3)



