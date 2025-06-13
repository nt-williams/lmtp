cf_sdr <- function(task, density_ratios, learners, control, progress_bar) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_sdr(task, fold, density_ratios, learners, control, progress_bar)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
       shifted = recombine(rbind_depth(ans, "shifted"), task$folds),
       fits = lapply(ans, function(x) x[["fits"]]))
}

estimate_sdr <- function(task, fold, density_ratios, learners, control, progress_bar) {
  # Get data splits for the current fold
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  density_ratios <- get_folded_data(density_ratios, task$folds, fold)$train

  # Pre-allocate matrices to store predictions
  pred_natural_train <- matrix(nrow = nrow(natural$train), ncol = task$time_horizon + 1)
  pred_shifted_train <- pred_natural_train

  pred_natural_valid <- matrix(nrow = nrow(natural$valid), ncol = task$time_horizon + 1)
  pred_shifted_valid <- pred_natural_valid

  pred_shifted_train[, task$time_horizon + 1] <- natural$train[[task$vars$Y]]
  pred_shifted_valid[, task$time_horizon + 1] <- natural$valid[[task$vars$Y]]

  # Pre-allocate list to store fitted models
  fits <- vector("list", length = task$time_horizon)

  # Loop over time points in reverse order
  for (time in rev(seq_len(task$time_horizon))) {
    y1 <- task$is_outcome_free(natural$train, time - 1)
    d0 <- task$is_competing_risk_free(natural$train, time - 1)
    c1 <- task$observed(natural$train, time)
    i <- c1 %and% (y1 & d0)

    history <- task$vars$history("L", time + 1)
    vars <- c("..i..lmtp_id", history, task$vars$Y)

    fit <- run_ensemble(natural$train[i, vars], task$vars$Y,
                        learners,
                        ifelse(time != task$time_horizon, "continuous", task$outcome_type),
                        "..i..lmtp_id",
                        control$.learners_outcome_folds)

    if (control$.return_full_fits) {
      fits[[time]] <- fit
    } else {
      fits[[time]] <- extract_sl_weights(fit)
    }

    A_t <- current_trt(task$vars$A, time)

    cp1 <- task$observed(natural$train, time - 1) # censoring in the past = 1
    y1v <- task$is_outcome_free(natural$valid, time - 1)
    d0v <- task$is_competing_risk_free(natural$valid, time - 1)
    cp1v <- task$observed(natural$valid, time - 1)

    i <- cp1 %and% (y1 & d0)
    iv <- cp1v %and% (y1v & d0v)

    under_shift_train <- natural$train[i, c("..i..lmtp_id", history)]
    under_shift_train[, A_t] <- shifted$train[i, A_t]

    pred_natural_train[i, time] <- predict(fit, natural$train[i, ], NULL)
    pred_shifted_train[i, time] <- predict(fit, under_shift_train, NULL)

    pred_natural_train[which(!y1), time] <- 0
    pred_natural_train[which(!d0), time] <- 1
    pred_shifted_train[which(!y1), time] <- 0
    pred_shifted_train[which(!d0), time] <- 1

    under_shift_valid <- natural$valid[iv, c("..i..lmtp_id", history)]
    under_shift_valid[, A_t] <- shifted$valid[iv, A_t]

    pred_natural_valid[iv, time] <- predict(fit, natural$valid[iv, ], NULL)
    pred_shifted_valid[iv, time] <- predict(fit, under_shift_valid, NULL)

    pred_natural_valid[which(!y1v), time] <- 0
    pred_natural_valid[which(!d0v), time] <- 1
    pred_shifted_valid[which(!y1v), time] <- 0
    pred_shifted_valid[which(!d0v), time] <- 1

    natural$train[, task$vars$Y] <- eif(density_ratios, pred_shifted_train, pred_natural_train, time)

    progress_bar()
  }

  list(natural = pred_natural_valid,
       shifted = pred_shifted_valid,
       fits = fits)
}
