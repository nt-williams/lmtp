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
  estimates <- fix_surv_time1(estimates)
  estimates <- isotonic_projection(estimates)

  class(estimates) <- "lmtp_survival"
  estimates
}

isotonic_projection <- function(x, alpha = 0.05) {
  cv <- abs(qnorm(p = alpha / 2))
  estim <- do.call("rbind", lapply(x, tidy))
  iso_fit <- isotone::gpava(1:length(x), estim$estimate)
  for (i in seq_along(x)) {
    x[[i]]$theta <- iso_fit$y[i]
    x[[i]]$low  <- x[[i]]$theta - (qnorm(0.975) * x[[i]]$standard_error)
    x[[i]]$high <- x[[i]]$theta + (qnorm(0.975) * x[[i]]$standard_error)
  }
  x
}