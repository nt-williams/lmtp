estimate_tmle_delay <- function(task, fold, propensity_scores, learners, control, progress_bar) {
  data <- get_folded_data(task$natural, task$folds, fold)
  levels <- task$support(task$time_horizon)
  number_levels <- length(levels)

  propensity_scores <- get_folded_data(propensity_scores, task$folds, fold)$train

  # Create empty array to store predictions
  # TODO: better logic column names, currently assumes the support is the same over time
  predictions_train_m <- array(NA_real_, c(nrow(data$train), number_levels, task$time_horizon), dimnames = list(NULL, levels, NULL))
  predictions_train_q <- predictions_train_m

  predictions_valid_m <- array(NA_real_, c(nrow(data$valid), number_levels, task$time_horizon), dimnames = list(NULL, levels, NULL))
  predictions_valid_q <- predictions_valid_m

  # Loop backwards in time for sequential regressions
  for (time in rev(seq_len(task$time_horizon))) {
    browser()
    y1 <- task$is_outcome_free(data$train, time - 1)
    d0 <- task$is_competing_risk_free(data$train, time - 1)
    c1 <- task$observed(data$train, time)
    i <- c1 %and% (y1 & d0)

    history <- task$vars$history("L", time + 1)
    vars <- c("..i..lmtp_id", history, task$vars$Y)

    # Estimate the outcome regressions
    # If its the last time point just perform this once using the real outcome
    if (time == task$time_horizon) {
      fit <- run_ensemble(
        data$train[i, vars], task$vars$Y, learners,
        task$outcome_type, "..i..lmtp_id", control$.learners_outcome_folds
      )

      learner_summary <- summary(fit, time, fold, level = NA_character_)
    }

    if (time < task$time_horizon) {
      # Need to fit regressions for each level of support
      # Create list to store fits for each level
      fits <- vector("list", length(task$support(time)))
      names(fits) <- task$support(time)

      for (s in task$support(time)) {
        data$train[[task$vars$Y]] <- predictions_train_q[, as.character(s), time + 1]

        fits[[as.character(s)]] <- run_ensemble(
          data$train[i, vars], task$vars$Y, learners,
          "continuous", "..i..lmtp_id", control$.learners_outcome_folds
        )

        learner_summary <- rbind(learner_summary, summary(fits[[as.character(s)]], time, fold, s))
      }
    }

    # Get treatment at this time
    this_treatment <- current_trt(task$vars$A, time)

    # Create subset indicators for survival, competing risk, censoring
    cp1 <- task$observed(data$train, time - 1)
    y1v <- task$is_outcome_free(data$valid, time - 1)
    d0v <- task$is_competing_risk_free(data$valid, time - 1)
    cp1v <- task$observed(data$valid, time - 1)

    ip <- cp1 %and% (y1 & d0)
    iv <- cp1v %and% (y1v & d0v)

    # Generate predictions
    if (time == task$time_horizon) {
      for (s in task$support(time - 1)) {
        # Predictions for validation data
        tmp <- data$valid
        tmp[[this_treatment]] <- s
        predictions_valid_q[iv, as.character(s), time] <- predict(fit, tmp[iv, ], 1e-05)

        # Predictions for training data
        tmp <- data$train
        tmp[[this_treatment]] <- s
        predictions_train_q[ip, as.character(s), time] <- predict(fit, tmp[ip, ], 1e-05)
      }
    }

    if (time < task$time_horizon) {
      # Loop over values to set
      for (s in task$support(time - 1)) {
        # Predictions for training
        tmp <- data$train
        tmp[[this_treatment]] <- s
        this_treatment_vector <- data$train[ip, this_treatment]
        pred <- vector("numeric", length(this_treatment_vector))

        # Loop over fits
        for (s2 in task$support(time)) {
          pred[this_treatment_vector == s2] <- predict(fits[[as.character(s2)]], tmp[ip, ], 1e-05)[this_treatment_vector == s2]
        }
        predictions_train_q[ip, as.character(s), time] <- pred

        # Predictions for validation
        tmp <- data$valid
        tmp[[this_treatment]] <- s
        this_treatment_vector <- data$valid[iv, this_treatment]
        pred <- vector("numeric", length(this_treatment_vector))

        for (s2 in task$support(time)) {
          pred[this_treatment_vector == s2] <- predict(fits[[as.character(s2)]], tmp[iv, ], 1e-05)[this_treatment_vector == s2]
        }
        predictions_valid_q[iv, as.character(s), time] <- pred
      }
    }

    # Handle deterministic predictions from survival/competing risks
    predictions_train_q[which(!y1), , time] <- 0
    predictions_train_q[which(!d0), , time] <- 1
    predictions_valid_q[which(!y1v), , time] <- 0
    predictions_valid_q[which(!d0v), , time] <- 1

    # TODO: iterate the progress bar
  }

  predictions_valid_q
}
