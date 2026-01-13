estimate_tmle_delay <- function(task, fold, propensity_scores, learners, control, progress_bar) {
  data <- get_folded_data(task$natural, task$folds, fold)

  # Augment data
  train <- delay_augment(data$train, task$sequences(task$time_horizon))
  valid <- delay_augment(data$valid, task$sequences(task$time_horizon))

  propensity_scores <- get_folded_data(propensity_scores, task$folds, fold)$train

  # Pre-allocate list to store predictions
  predictions <- vector("list", task$time_horizon)

  # Loop backwards in time for sequential regressions
  for (time in rev(seq_len(task$time_horizon))) {
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

    # Treatment sequences we need to generate under
    sequences <- task$sequences(time)

    if (time < task$time_horizon) {
      # Fit regressions for each level of support using a pooled regression
      # Add the pooling variable
      vars <- c(vars, paste0("..i..lmtp_tmp_s", time))
      fit <- run_ensemble(
        subset_augmented(train, time, task$time_horizon)[delay_augment(i, sequences), vars],
        task$vars$Y, learners,
        "continuous", "..i..lmtp_id",
        control$.learners_outcome_folds
      )

      learner_summary <- rbind(learner_summary, summary(fit, time, fold, level = NA_character_))
    }

    # Generate predictions
    # Get treatment at this time
    this_treatment <- current_trt(task$vars$A, time)

    # Create subset indicators for survival, competing risk, censoring
    cp1 <- task$observed(data$train, time - 1)
    y1v <- task$is_outcome_free(data$valid, time - 1)
    d0v <- task$is_competing_risk_free(data$valid, time - 1)
    cp1v <- task$observed(data$valid, time - 1)

    ip <- cp1 %and% (y1 & d0)
    iv <- cp1v %and% (y1v & d0v)

    # Store Y_(t+1) before overriding
    outcome <- subset_augmented(train, time, task$time_horizon)[[task$vars$Y]]

    # Generate predictions
    train <- predict_delay_augment(
      subset_augmented(train, time, task$time_horizon), fit,
      sequences, time, task$time_horizon,
      this_treatment, task$vars$Y,
      ip, y1, d0
    )

    valid <- predict_delay_augment(
      subset_augmented(valid, time, task$time_horizon), fit,
      sequences, time, task$time_horizon,
      this_treatment, task$vars$Y,
      iv, y1v, d0v
    )

    weights <- create_tmle_delay_weights(
      train,
      task$vars$A, this_treatment,
      time,
      propensity_scores
    )

    fit <- fluc(outcome[delay_augment(i, sequences)],
                train[delay_augment(i, sequences), task$vars$Y],
                weights[delay_augment(i, sequences)])

    train[delay_augment(i, sequences), task$vars$Y] <-
      update(fit, train[delay_augment(i, sequences), task$vars$Y])
    valid[delay_augment(iv, sequences), task$vars$Y] <-
      update(fit, valid[delay_augment(iv, sequences), task$vars$Y])

    # Save all validation predictions
    predictions[[time]] <- valid[, c("..i..lmtp_id", paste0("..i..lmtp_tmp_s", seq_len(time)), task$vars$Y)]

    # TODO: iterate the progress bar
  }

  predictions
}
