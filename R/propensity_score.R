cf_propensity_score <- function(task, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_propensity_score(task, fold, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  ans <- list(
    propensity_score = lapply(seq_len(task$time_horizon), function(time) {
      this_time <- lapply(seq_along(task$folds), function(fold) {
        ans[[fold]]$propensity_score[[time]]
      })
      recombine(Reduce(rbind, this_time), task$folds)
    }),
    prob_observed = recombine(
      rbind_depth(ans, "prob_observed"),
      task$folds
    ),
    fits_trt = lapply(ans, function(x) x[["fits_trt"]]),
    fits_cens = lapply(ans, function(x) x[["fits_cens"]])
  )

  # TODO: trimming
  # ans$propensity_score <- trim(ans$propensity_score, control$.trim)
  ans
}

estimate_propensity_score <- function(task, fold, learners, control, pb) {
  data <- get_folded_data(task$natural, task$folds, fold)

  # Create list of matrices the same number as \tau to store predictions
  propensity_scores <- lapply(seq_len(task$time_horizon),
                              \(time) matrix(nrow = nrow(data$valid), ncol = 2))
  prob_observed <- matrix(1, nrow = nrow(data$valid), ncol = task$time_horizon)
  fits_trt <- vector("list", length = task$time_horizon)
  fits_cens <- fits_trt

  for (time in seq_len(task$time_horizon)) {
    # Get indices of observations that aren't censored and haven't experienced
    # the outcome or the competing risk
    i <- task$observed(data$train, time - 1) %and%
      task$is_at_risk(data$train, time)

    # Current treatment variable
    A_t <- current_trt(task$vars$A, time)
    # Current censoring variable
    C_t <- task$vars$C[time]

    # Covariates
    vars <- c("..i..lmtp_id", task$vars$history("A", time))

    # Need to fit K-1 models for K levels
    # Just assuming binary for now
    fit_trt <- run_ensemble(
      data$train[i, c(vars, A_t)], A_t,
      learners,
      "binomial", "..i..lmtp_id",
      control$.learners_trt_folds
    )

    if (!is.null(C_t)) {
      fit_cens <- run_ensemble(
        data[i, c(vars, A_t, C_t)], C_t,
        learners,
        "binomial", "..i..lmtp_id",
        control$.learners_trt_folds
      )
    } else {
      fit_cens <- NULL
    }

    if (control$.return_full_fits) {
      fits_trt[[time]] <- fit_trt
      fits_cens[[time]] <- fit_cens
    } else {
      fits_trt[[time]] <- extract_sl_weights(fit_trt)
      fits_cens[[time]] <- extract_sl_weights(fit_cens)
    }

    i <- task$observed(data$valid, time - 1) %and%
      task$is_at_risk(data$valid, time)

    # j=2 will index A = 1, j=1 will index A = 0
    # TODO: dont hard code this
    propensity_scores[[time]][i, 2] <- predict(fit_trt, data$valid[i, ])
    propensity_scores[[time]][i, 1] <- 1 - propensity_scores[[time]][i, 2]
    colnames(propensity_scores[[time]]) <- c("0", "1")

    if (!is.null(C_t)) {
      # TODO: predict probability of being observed
    }

    pb()
  }

  list(propensity_score = propensity_scores,
       prob_observed = prob_observed,
       fits_trt = fits_trt,
       fits_cens = fits_cens)
}
