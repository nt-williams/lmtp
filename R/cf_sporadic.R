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
  # create folds from the "natural" data
  natural <- get_folded_data(task$natural, task$folds, fold)

  # initialize matrix to store predictions, one column for each time point initially set to 1
  weights <- matrix(data = 1, nrow = nrow(natural$valid), ncol = task$tau)
  # list to store model fits
  fits <- vector("list", length = task$tau)

  # for now this only programmed to work with a survival outcome
  for (t in 1:(task$tau - 1)) {
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
    ..i..R <- as.numeric(natural$train[[task$vars$C[t]]] == 1 & is.na(natural$train[[task$vars$N[t]]]))

    # Adding missingess indicator to training data
    natural$train$..i..R <- ..i..R

    # Variables needed for estimation
    vars <- c("..i..lmtp_id", task$vars$history("A", t), A_t, "..i..R")

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

    # save fit
    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    # Indicators for censoring in the validation set
    i <- task$observed(natural$valid, t - 1)

    # Create an indicator for sporadic missingness in the validation set
    ..i..R <- natural$valid[[task$vars$C[t]]] == 1 & is.na(natural$valid[[task$vars$N[t]]])

    # Predict on validation set
    pred <- matrix(1, nrow = nrow(natural$valid), ncol = 1)
    pred[i, ] <- 1 - predict(fit, natural$valid[i, ])

    # Create IPW weights
    weights[, t] <- (!..i..R) / pred
  }

  list(weights = weights, fits = fits)
}
