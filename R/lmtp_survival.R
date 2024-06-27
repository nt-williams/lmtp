lmtp_survival <- function(data, trt, outcomes, baseline = NULL, time_vary = NULL,
  cens = NULL, shift = NULL, shifted = NULL, estimator = c("lmtp_tmle", "lmtp_sdr"), k = Inf,
  mtp = FALSE, 
  id = NULL, 
  learners_outcome = "SL.glm",
  learners_trt = "SL.glm",
  folds = 10, weights = NULL,
  control = lmtp_control()) {
    
  checkmate::assertCharacter(outcomes, min.len = 2, null.ok = FALSE, unique = TRUE, any.missing = FALSE)

  estimator <- match.arg(estimator)
  tau <- length(outcomes)
  estimates <- vector("list", tau)

  for (t in 1:tau) {
    args <- list(
      data = data, 
      trt = trt, 
      outcome = outcomes[1:t], 
      baseline = baseline, 
      time_vary = time_vary, 
      cens = cens[1:t], 
      shift = shift, 
      shifted = shifted, 
      k = k, 
      mtp = mtp, 
      outcome_type = ifelse(t == 1, "binomial", "survival"), 
      id = id, 
      learners_outcome = learners_outcome, 
      learners_trt = learners_trt, 
      folds = folds, 
      weights = weights, 
      control = control
    )

    estimates[[t]] <- future::future({
      if (estimator == "lmtp_tmle") do.call(lmtp_tmle, args)
      else do.call(lmtp_sdr, args)
    },
    seed = TRUE)
  }

  estimates <- future::value(estimates)

  estimates[[1]]$theta <- 1 - estimates[[1]]$theta
  high <- estimates[[1]]$high
  low <- estimates[[1]]$low
  estimates[[1]]$high <- 1 - low
  estimates[[1]]$low <- 1 - high
  estimates[[1]]$eif <- 1 - estimates[[1]]$eif

  class(estimates) <- "lmtp_survival"
  estimates
}