cf_sporadic <- function(task, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_sporadic(task, fold, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(weights = recombine(rbind_depth(ans, "weights"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_sporadic <- function(task, fold, learners, control, pb) {
  # Create folds from the "natural" data
  natural <- get_folded_data(task$natural, task$folds, fold)

  # Initialize matrix to store predictions, one column for each time point initially set to 1
  weights <- matrix(data = 1, nrow = nrow(natural$valid), ncol = task$tau)
  # List to store model fits
  fits <- vector("list", length = task$tau)

  # For now this only programmed to work with a survival outcome
  # TODO: NEED TO START PASSING TIME-VARYING OUTCOMES TO THE OUTCOME NODE
  for (t in 1:(task$tau - 1)) {
    browser()
    # Indicators for censoring
    i <- task$observed(natural$train, t - 1)

    # Figure out the treatment at time t
    if (length(task$vars$A) > 1) {
      # if there are multiple treatments, we need to use the treatment at time t
      A_t <- task$vars$A[[t]]
    } else {
      # if there is only one treatment, we can use the treatment at time 1
      A_t <- task$vars$A[[1]]
    }

    # Create an indicator for sporadic missingness
    # If there is no time-varying outcome, there can't be sporadic outcome measurement
    if (is.null(task$vars$N[t])) {
      ..i..R <- rep(0, nrow(natural$train))
    } else {
      # Otherwise, sporadic outcome measurement exists if an observation is uncensored, but the outcome is missing
      ..i..R <- as.numeric(task$observed(natural$train, t) & is.na(natural$train[[task$vars$N[t]]]))
    }

    # Adding sporadic measurement indicator to training data
    natural$train$..i..R <- ..i..R

    # Variables needed for estimation
    vars <- c("..i..lmtp_id", task$vars$history("A", t), A_t, "..i..R")

    # Simple test that there is variation in the sporadic measurement indicator
    if (var(..i..R) < 0.0001) {
      fit <- lm(..i..R ~ 1, data = natural$train[i, vars])
    } else {
      # Fit model for sporadic missingness using the super learner
      fit <- run_ensemble(natural$train[i, vars], "..i..R",
                          learners, "binomial", "..i..lmtp_id",
                          control$.learners_trt_folds,
                          control$.discrete,
                          control$.info)
    }

    # Save fit
    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    # Indicators for censoring in the validation set
    i <- task$observed(natural$valid, t - 1)

    # Create an indicator for sporadic measurement in the validation set
    # Using same logic as for training data
    if (is.null(task$vars$N[t])) {
      ..i..R <- rep(0, nrow(natural$train))
    } else {
      ..i..R <- as.numeric(task$observed(natural$valid, t) & is.na(natural$valid[[task$vars$N[t]]]))
    }

    # Predict on validation set
    pred <- matrix(1, nrow = nrow(natural$valid), ncol = 1)
    pred[i, ] <- 1 - predict(fit, natural$valid[i, ])

    # Create IPW weights
    weights[, t] <- (!..i..R) / pred
  }

  list(weights = weights, fits = fits)
}
