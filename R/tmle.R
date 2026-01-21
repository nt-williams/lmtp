cf_tmle <- function(task, density_ratios, learners, control, progress_bar) {
  ans <- vector("list", length = length(task$folds))

  density_ratios <- compute_weights(density_ratios, 1, task$time_horizon)

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_tmle(task, fold, density_ratios, learners, control, progress_bar)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(predictions = recombine(rbind_depth(ans, "predictions"), task$folds),
       uncentered_eif = recombine(c_depth(ans, "uncentered_eif"), task$folds),
       learner_outcome_summary = rbind_depth(ans, "learner_summary"))
}

estimate_tmle <- function(task, fold, density_ratios, learners, control, progress_bar) {
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

  fits <- vector("list", length = time_horizon)
  for (time in rev(seq_len(time_horizon))) {
    y1 <- task$is_outcome_free(natural_train, time - 1)
    d0 <- task$is_competing_risk_free(natural_train, time - 1)
    c1 <- task$observed(natural_train, time)
    i <- c1 %and% (y1 & d0)

    history <- task$vars$history("L", time + 1)
    vars <- c(id, history, Y)

    fit <- run_ensemble(
      natural_train[i, vars], Y, learners,
      ifelse(time != task$time_horizon, "continuous", task$outcome_type),
      id, cvfolds
    )

    learner_summaries[[time]] <- summary(fit, time, fold)

    A_t <- current_trt(A, time)

    cp1 <- task$observed(natural_train, time - 1)
    y1v <- task$is_outcome_free(natural_valid, time - 1)
    d0v <- task$is_competing_risk_free(natural_valid, time - 1)
    cp1v <- task$observed(natural_valid, time - 1)

    ip <- cp1 %and% (y1 & d0)
    iv <- cp1v %and% (y1v & d0v)

    under_shift_train <- natural_train[ip, c(id, history)]
    under_shift_train[, A_t] <- shifted_train[ip, A_t]

    pred_natural_train[ip, time] <- predict(fit, natural_train[ip, ], 1e-05)
    pred_shifted_train[ip, time] <- predict(fit, under_shift_train, 1e-05)

    under_shift_valid <- natural_valid[iv, c(id, history)]
    under_shift_valid[, A_t] <- shifted$valid[iv, A_t]

    pred_natural_valid[iv, time] <- predict(fit, natural_valid[iv, ], 1e-05)
    pred_shifted_valid[iv, time] <- predict(fit, under_shift_valid, 1e-05)

    # Fit fluctuation model and update predictions
    fit <- fluc(natural_train[i, Y],
                pred_natural_train[i, time],
                densrat_train[i, time] * weights$train[i])

    natural_train[ip, Y] <- update(fit, pred_shifted_train[ip, time])

    pred_natural_valid[iv, time] <- update(fit, pred_natural_valid[iv, time])
    pred_shifted_valid[iv, time] <- update(fit, pred_shifted_valid[iv, time])

    # Calibrate deterministic predictions
    natural_train[[Y]] <- calibrate(natural_train[[Y]], y1, d0)
    pred_natural_valid[, time] <- calibrate(pred_natural_valid[, time], y1v, d0v)
    pred_shifted_valid[, time] <- calibrate(pred_shifted_valid[, time], y1v, d0v)

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

fluc <- function(y, offset, weights) {
  sw(glm(y ~ offset(qlogis(offset)), weights = weights, family = "binomial"))
}

update <- function(fluc, pred) {
  bound(plogis(qlogis(pred) + coef(fluc)))
}
