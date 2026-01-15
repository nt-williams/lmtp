cf_tmle_delay <- function(task, propensity_scores, learners, control, progress_bar) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_tmle_delay(task, fold, propensity_scores, learners, control, progress_bar)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(predictions = recombine(c_depth(ans, "predictions"), task$folds),
       uncentered_eif = recombine(c_depth(ans, "uncentered_eif"), task$folds),
       learner_outcome_summary = rbind_depth(ans, "learner_summary"))
}

estimate_tmle_delay <- function(task, fold, propensity_scores, learners, control, progress_bar) {
  data <- get_folded_data(task$natural, task$folds, fold)

  # Cache used objects to prevent look-ups
  time_horizon <- task$time_horizon
  seq_time_horizon <- task$sequences(time_horizon)

  # Augment data
  train <- data$train
  valid <- data$valid
  train_aug <- delay_augment(train, seq_time_horizon)
  valid_aug <- delay_augment(valid, seq_time_horizon)

  propensity_scores <- get_folded_data(propensity_scores, task$folds, fold)

  # Pre-allocate lists
  pred_train_shifted <- vector("list", time_horizon)
  pred_valid_shifted <- pred_train_shifted
  learner_summaries <- vector("list", time_horizon)

  # Pre-allocate vector for the EIF
  eif <- numeric(nrow(valid))

  # Cache used objects to prevent look-ups
  outcomevar <- task$vars$Y
  trtvars <- task$vars$A
  idvar <- "..i..lmtp_id"
  cvfolds <- control$.learners_outcome_folds

  # Pre-computations
  all_sequences <- lapply(seq_len(time_horizon), task$sequences)

  # Loop backwards in time for sequential regressions
  for (time in rev(seq_len(time_horizon))) {
    y1 <- task$is_outcome_free(train, time - 1)
    d0 <- task$is_competing_risk_free(train, time - 1)
    c1 <- task$observed(train, time)
    i <- c1 %and% (y1 & d0)

    # Cache used objects to prevent look-ups
    history <- task$vars$history("L", time + 1)
    vars <- c(idvar, history, outcomevar)

    # Estimate the outcome regressions
    # If its the last time point just perform this once using the real outcome
    if (time == task$time_horizon) {
      fit <- run_ensemble(train[i, vars], outcomevar, learners, task$outcome_type, idvar, cvfolds)
      learner_summaries[[time]] <- summary(fit, time, fold)
    }

    # Treatment sequences we need to generate under
    sequences <- all_sequences[[time]]
    i_aug <- delay_augment(i, sequences)

    if (time < time_horizon) {
      train_aug[[outcomevar]] <- pred_train_shifted[[time + 1]]
      valid_aug[[outcomevar]] <- pred_valid_shifted[[time + 1]]
    }

    # Subset the augmented data
    train_aug <- subset_augmented(train_aug, time, time_horizon)
    valid_aug <- subset_augmented(valid_aug, time, time_horizon)

    if (time < time_horizon) {
      # Fit regressions for each level of support using a pooled regression
      # Add the pooling variable
      vars <- c(vars, seqvars(time))
      fit <- run_ensemble(train_aug[i_aug, vars], outcomevar, learners, "continuous", idvar, cvfolds)
      learner_summaries[[time]] <- summary(fit, time, fold)
    }

    # Generate predictions
    # Get treatment at this time
    this_treatment <- current_trt(trtvars, time)

    # Create subset indicators for survival, competing risk, censoring
    cp1 <- task$observed(train, time - 1)
    y1v <- task$is_outcome_free(valid, time - 1)
    d0v <- task$is_competing_risk_free(valid, time - 1)
    cp1v <- task$observed(valid, time - 1)

    ip <- cp1 %and% (y1 & d0)
    ip_aug <- delay_augment(ip, sequences)
    y1_aug <- delay_augment(y1, sequences)
    d0_aug <- delay_augment(d0, sequences)

    iv <- cp1v %and% (y1v & d0v)
    iv_aug <- delay_augment(iv, sequences)
    y1v_aug <- delay_augment(y1v, sequences)
    d0v_aug <- delay_augment(d0v, sequences)

    # Generate predictions
    pred_train_shifted[[time]] <- predict_delay_augment(
      train_aug, fit, time, time_horizon, this_treatment,
      outcomevar, ip_aug, y1_aug, d0_aug, TRUE
    )

    pred_train_natural <- predict_delay_augment(
      train_aug, fit, time, time_horizon, this_treatment,
      outcomevar, ip_aug, y1_aug, d0_aug, FALSE
    )

    pred_valid_shifted[[time]] <- predict_delay_augment(
      valid_aug, fit, time, time_horizon, this_treatment,
      outcomevar, iv_aug, y1v_aug, d0v_aug, TRUE
    )

    pred_valid_natural <- predict_delay_augment(
      valid_aug, fit, time, time_horizon, this_treatment,
      outcomevar, iv_aug, y1v_aug, d0v_aug, FALSE
    )

    # Riesz representer
    # weights <- tmle_delay_weights(train_aug, trtvars, this_treatment, time, propensity_scores$train)
    weights <- delay_riesz_rep(train_aug, trtvars, propensity_scores$train, this_treatment, 1, time)

    # Fit tilting model
    fit <- fluc(train_aug[[outcomevar]][i_aug], pred_train_natural[i_aug], weights[i_aug])

    # Update predictions
    pred_train_shifted[[time]][i_aug] <- update(fit, pred_train_shifted[[time]][i_aug])
    pred_valid_shifted[[time]][iv_aug] <- update(fit, pred_valid_shifted[[time]][iv_aug])
    pred_valid_natural[iv_aug] <- update(fit, pred_valid_natural[iv_aug])

    # Construct the EIF
    weights <- tmle_delay_weights(valid_aug, trtvars, this_treatment, time, propensity_scores$valid)
    eif_comp <- weights * (valid_aug[[outcomevar]] - pred_valid_natural)
    eif <- eif + collapse::fsum(eif_comp, valid_aug[[idvar]])

    # TODO: iterate the progress bar
  }

  # Time 0 component
  valid_aug[[outcomevar]] <- pred_valid_shifted[[1]]
  valid_aug <- subset_augmented(valid_aug, 0, time_horizon)
  # eif <- as.vector(eif + (valid_aug[[outcomevar]] - collapse::fmean(valid_aug[[outcomevar]])))
  eif <- as.vector(eif + valid_aug[[outcomevar]])

  list(predictions = valid_aug[[outcomevar]],
       uncentered_eif = eif,
       learner_summary = data.table::rbindlist(learner_summaries))
}
