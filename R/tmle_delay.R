estimate_tmle_delay <- function(task, fold, propensity_scores, learners, control, progress_bar) {
  data <- get_folded_data(task$natural, task$folds, fold)

  # Augment data
  train <- delay_augment(data$train, task$sequences(task$time_horizon))
  valid <- delay_augment(data$valid, task$sequences(task$time_horizon))

  propensity_scores <- get_folded_data(propensity_scores, task$folds, fold)$train

  # Pre-allocate list to store predictions
  pred_train_shifted <- vector("list", task$time_horizon)
  pred_valid_shifted <- pred_train_shifted

  ic <- numeric(nrow(data$valid))

  # Loop backwards in time for sequential regressions
  for (time in rev(seq_len(task$time_horizon))) {
    y1 <- task$is_outcome_free(data$train, time - 1)
    d0 <- task$is_competing_risk_free(data$train, time - 1)
    c1 <- task$observed(data$train, time)
    i <- c1 %and% (y1 & d0)

    history <- task$vars$history("L", time + 1)
    idvar <- "..i..lmtp_id"
    seqvar <- "..i..lmtp_tmp_s"
    vars <- c(idvar, history, task$vars$Y)

    # Estimate the outcome regressions
    # If its the last time point just perform this once using the real outcome
    if (time == task$time_horizon) {
      fit <- run_ensemble(
        data$train[i, vars], task$vars$Y, learners,
        task$outcome_type, idvar, control$.learners_outcome_folds
      )

      learner_summary <- summary(fit, time, fold, level = NA_character_)
    }

    # Treatment sequences we need to generate under
    sequences <- task$sequences(time)

    if (time < task$time_horizon) {
      train[[task$vars$Y]] <- pred_train_shifted[[time + 1]]
    }

    # Subset the augmented data
    train <- subset_augmented(train, time, task$time_horizon)
    valid <- subset_augmented(valid, time, task$time_horizon)

    if (time < task$time_horizon) {
      # Fit regressions for each level of support using a pooled regression
      # Add the pooling variable
      vars <- c(vars, paste0(seqvar, time))
      fit <- run_ensemble(
        train[delay_augment(i, sequences), vars],
        task$vars$Y, learners,
        "continuous", idvar,
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

    # Generate predictions
    pred_train_shifted[[time]] <- predict_delay_augment(
      train, fit, sequences, time,
      task$time_horizon, this_treatment, task$vars$Y, ip, y1, d0, TRUE
    )

    pred_train_natural <- predict_delay_augment(
      train, fit, sequences, time,
      task$time_horizon, this_treatment, task$vars$Y, ip, y1, d0, FALSE
    )

    pred_valid_shifted[[time]] <- predict_delay_augment(
      valid, fit, sequences, time,
      task$time_horizon, this_treatment, task$vars$Y, iv, y1v, d0v, TRUE
    )

    pred_valid_natural <- predict_delay_augment(
      train, fit, sequences, time,
      task$time_horizon, this_treatment, task$vars$Y, iv, y1v, d0v, FALSE
    )

    # Riesz representer
    weights <- tmle_delay_weights(train, task$vars$A, this_treatment, time, propensity_scores)

    # Fit tilting model
    i_aug <- delay_augment(i, sequences)
    iv_aug <- delay_augment(iv, sequences)

    fit <- fluc(train[[task$vars$Y]][i_aug], pred_train_natural[i_aug], weights[i_aug])

    # Update predictions
    pred_train_shifted[[time]][i_aug] <- update(fit, pred_train_shifted[[time]][i_aug])
    pred_valid_shifted[[time]][iv_aug] <- update(fit, pred_valid_shifted[[time]][iv_aug])
    pred_valid_natural[iv_aug] <- update(fit, pred_valid_natural[iv_aug])

    # Construct the EIF
    ic_comp <- weights * (train[[task$vars$Y]] - pred_valid_natural)
    ic <- ic + collapse::fsum(ic_comp, valid[[idvar]])

    # TODO: iterate the progress bar
  }

  valid[[task$vars$Y]] <- pred_valid_shifted[[1]]
  valid <- subset_augmented(valid, 0, task$time_horizon)
  ic <- as.vector(ic + (valid[[task$vars$Y]] - collapse::fmean(valid[[task$vars$Y]])))

  list(predictions = pred_valid_shifted,
       efficient_influence_function = ic,
       learner_summary = learner_summary)
}
