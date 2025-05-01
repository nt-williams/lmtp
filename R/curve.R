cf_curve <- function(task, ratios, sporadic_weights, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_curve_sdr(task, fold, ratios, sporadic_weights, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(influence_functions = recombine(rbind_depth(ans, "influence_functions"), task$folds),
       shifted = lapply(1:task$tau, \(t) recombine(rbind_depth_2(ans, "shifted", t), task$folds)),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_curve_sdr <- function(task, fold, ratios, sporadic_weights, learners, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  ratios <- get_folded_data(ratios, task$folds, fold)
  sporadic_weights <- get_folded_data(sporadic_weights, task$folds, fold)

  # Matrix to store predictions under the intervention for the training data
  mst <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$train), ncol = tau + 1))
  # Matrix to store predictions under the observed for the training data
  mnt <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$train), ncol = tau))
  # Matrix to store predictions under the intervention for the validation data
  msv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau + 1))
  # Matrix to store predictions under the observed for the validation data
  mnv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau))

  # Matrix to store final influence functions
  influence_functions <- matrix(nrow = nrow(natural$valid), ncol = task$tau)

  # Pivot into long format
  natural <- lapply(natural, \(x) pivot(x, task$vars))
  shifted <- lapply(shifted, \(x) pivot(x, task$vars))

  # Keep lagged variables as natural values for predictions
  # Under shift training
  ust <- natural$train
  As <- grep("^..i..A", names(ust), value = TRUE)
  ust[, As] <- shifted$train[, As]

  # Under shift valid
  usv <- natural$valid
  usv[, As] <- shifted$valid[, As]

  taus <- order(unique(natural$train$time))
  for (t in seq_along(taus)) {
    mst[[t]][, t + 1] <- subset(natural$train, time == taus[t])$..i..Y_1
    msv[[t]][, t + 1] <- subset(natural$valid, time == taus[t])$..i..Y_1
  }

  fits <- vector("list", length = task$tau)
  for (t in 1:task$tau) {
    if (isTRUE(task$survival)) {
      # Indicator for not experiencing the outcome or the competing risk already
      at_risk <- natural$train$..i..N == 1 & natural$train$..i..D_1 == 0
      at_risk[is.na(at_risk)] <- FALSE
    } else {
      at_risk <- rep(TRUE, nrow(natural$train))
    }
    # Indicator for having the outcome measured (i.e., no sporadic measurement)
    is_measured <- natural$train$..i..R == 1
    # Indicator for not having been censored
    is_observed <- natural$train$..i..C_1 == 1
    # Combining indicators for the training subset
    train_subset <- at_risk & is_measured & is_observed

    # Create indicators for the subset of rows to use for estimation at this time
    # This uses slightly different logic than the paper to not rely on a variable 's'
    time <- as.numeric(natural$train$time) <= rev(1:task$tau)[t]

    # Remove variables we don't need for estimation
    vars <- setdiff(names(natural$train), c("..i..C_1", "..i..C_1_lag", "..i..wide_id", "..i..N", "..i..D_1", "..i..R"))
    # Don't need the time variable for the last iteration
    if (t == task$tau) {
      vars <- setdiff(vars, "time")
    }

    fit <- run_ensemble(
      natural$train[train_subset & time, vars],
      "..i..Y_1",
      learners,
      ifelse(t != 1, "continuous", task$outcome_type),
      "..i..lmtp_id",
      control$.learners_outcome_folds,
      control$.discrete,
      control$.info
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    # Update the matrices storing all predictions
    mnt <- update_m(mnt, fit, natural$train, t, task$tau)
    mst <- update_m(mst, fit, ust, t, task$tau)
    mnv <- update_m(mnv, fit, natural$valid, t, task$tau)
    msv <- update_m(msv, fit, usv, t, task$tau)

    # Only When t (l in the paper) = 1 do we apply the sporadic weights
    pseudo <- list()
    for (s in t:task$tau) {
      tau <- ncol(mnt[[s]])

      if (s == t) {
        # If the smallest time point in the sequence from t:tau use validation data to calculate EIF for theta at time t
        influence_functions[, t] <- eif(ratios$valid, sporadic_weights$valid, msv[[s]], mnv[[s]], t = 1, tau = tau)
      } else {
        # Otherwise calculate the DR transformation using training values as the pseudo outcome for the next iteration
        dr_transformation <- eif(ratios$train, sporadic_weights$train, mst[[s]], mnt[[s]], t = (s - t + 1), tau = tau)
        pseudo[[length(pseudo) + 1]] <- dr_transformation
      }
    }

    # Add the pseudo outcomes to the training data
    pseudo <- unlist(pseudo)
    natural$train[, "..i..Y_1"] <- c(pseudo, rep(NA_real_, nrow(natural$train) - length(pseudo)))

    # After the first iteration, R becomes 1 for all observations
    # Sub-setting is only based on censoring for iterations > 1
    if (t == 1) {
      natural$train$..i..R <- 1
    }

    pb()
  }

  list(
    influence_functions = influence_functions,
    shifted = msv,
    fits = fits,
    lmtp_id = natural$valid$..i..lmtp_id
  )
}

predict_long <- function(fit, newdata, t, tau) {
  # Create indicators for the subset of rows to use based on time
  time <- as.numeric(newdata$time) <= rev(1:tau)[t]
  # Empty matrix to store predictions
  predictions <- matrix(nrow = nrow(newdata[time, ]), ncol = 1)
  # Indicator for not having been censored at the previous time point
  is_observed <- newdata$..i..C_1_lag == 1

  if (isTRUE(task$survival)) {
    # Indicator for not experiencing the outcome already
    outcome_free <- newdata$..i..N[time] == 1
    # Indicator for not experiencing competing risk already
    competing_risk_free <- newdata$..i..D_1[time] == 0
  } else {
    outcome_free <- rep(TRUE, nrow(newdata[time, ]))
    competing_risk_free <- rep(TRUE, nrow(newdata[time, ]))
  }

  predictions[is_observed[time], 1] <- predict(fit, newdata[time & is_observed, ], NULL)
  predictions[!outcome_free, 1] <- 0
  predictions[!competing_risk_free, 1] <- 1
  predictions[, 1]
}

update_m <- function(m, fit, newdata, t, tau) {
  pred <- predict_long(fit, newdata, t, tau)
  time <- as.numeric(newdata$time) >= t

  for (l in seq_along(t:tau)) {
    x <- (t:tau)[l]
    j <- as.numeric(newdata$time[time]) == x
    m[[x]][, x - t + 1] <- pred[j]
  }
  m
}
