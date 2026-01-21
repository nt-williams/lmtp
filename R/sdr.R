cf_sdr <- function(task, density_ratios, learners, control, progress_bar) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_sdr(task, fold, density_ratios, learners, control, progress_bar)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(predictions = recombine(rbind_depth(ans, "predictions"), task$folds),
       uncentered_eif = recombine(c_depth(ans, "uncentered_eif"), task$folds),
       learner_outcome_summary = rbind_depth(ans, "learner_summary"))
}

estimate_sdr <- function(task, fold, density_ratios, learners, control, progress_bar) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  densrat <- get_folded_data(density_ratios, task$folds, fold)
  weights <- get_folded_data(task$weights, task$folds, fold)

  # Caching
  time_horizon <- task$time_horizon

  natural_train <- natural$train
  natural_valid <- natural$valid
  shifted_train <- shifted$train
  shifted_valid <- shifted$valid
  densrat_train <- densrat$train
  densrat_valid <- densrat$valid

  Y <- task$vars$Y
  A <- task$vars$A
  C <- task$vars$C
  id <- "..i..lmtp_id"
  cvfolds <- control$.learners_outcome_folds

  learner_summaries <- vector("list", time_horizon)

  pred_natural_train <- matrix(nrow = nrow(natural_train),
                               ncol = task$time_horizon + 1)
  pred_shifted_train <- matrix(nrow = nrow(shifted_train),
                               ncol = task$time_horizon + 1)

  pred_natural_valid <- matrix(nrow = nrow(natural_valid),
                               ncol = task$time_horizon + 1)
  pred_shifted_valid <- matrix(nrow = nrow(shifted_valid),
                               ncol = task$time_horizon + 1)

  pred_shifted_valid[, time_horizon + 1] <- natural_valid[[Y]]

  # Pre-allocate vector for the EIF
  eif <- numeric(nrow(natural_valid))
  # Cumprod valid density ratios for the EIF
  densrat_valid <- compute_weights(densrat_valid, 1, time_horizon)

  # Loop over time points in reverse order
  for (time in rev(seq_len(task$time_horizon))) {
    y1 <- task$is_outcome_free(natural_train, time - 1)
    d0 <- task$is_competing_risk_free(natural_train, time - 1)
    c1 <- task$observed(natural_train, time)
    i <- c1 %and% (y1 & d0)

    history <- task$vars$history("L", time + 1)
    vars <- c(id, history, Y)

    fit <- run_ensemble(
      natural_train[i, vars], Y, learners,
      ifelse(time != time_horizon, "continuous", task$outcome_type),
      id, cvfolds
    )

    learner_summaries[[time]] <- summary(fit, time, fold)

    A_t <- current_trt(A, time)

    cp1 <- task$observed(natural_train, time - 1) # censoring in the past = 1
    y1v <- task$is_outcome_free(natural_valid, time - 1)
    d0v <- task$is_competing_risk_free(natural_valid, time - 1)
    cp1v <- task$observed(natural_valid, time - 1)

    i <- cp1 %and% (y1 & d0)
    iv <- cp1v %and% (y1v & d0v)

    under_shift_train <- natural_train[i, c(id, history)]
    under_shift_train[, A_t] <- shifted_train[i, A_t]

    pred_natural_train[i, time] <- predict(fit, natural_train[i, ], NULL)
    pred_shifted_train[i, time] <- predict(fit, under_shift_train, NULL)

    # Calibrate deterministic predictions
    pred_natural_train[, time] <- calibrate(pred_natural_train[, time], y1, d0)
    pred_shifted_train[, time] <- calibrate(pred_shifted_train[, time], y1, d0)

    under_shift_valid <- natural_valid[iv, c(id, history)]
    under_shift_valid[, A_t] <- shifted_valid[iv, A_t]

    pred_natural_valid[iv, time] <- predict(fit, natural_valid[iv, ], NULL)
    pred_shifted_valid[iv, time] <- predict(fit, under_shift_valid, NULL)

    # Calibrate deterministic predictions
    pred_natural_valid[, time] <- calibrate(pred_natural_valid[, time], y1v, d0v)
    pred_shifted_valid[, time] <- calibrate(pred_shifted_valid[, time], y1v, d0v)

    # Update outcome with SDR transformation
    natural_train[, Y] <- eif(densrat_train, pred_shifted_train, pred_natural_train, time)

    # Construct the EIF
    eif <- update_lmtp_eif(
      eif, pred_shifted_valid[, time + 1],
      densrat_valid[, time], pred_natural_valid[, time]
    )

    progress_bar()
  }

  # EIF time 0 component
  eif <- eif + pred_shifted_valid[, 1]

  list(predictions = pred_shifted_valid,
       uncentered_eif = eif,
       learner_summary = data.table::rbindlist(learner_summaries))
}
