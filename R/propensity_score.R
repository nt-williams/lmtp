cf_propensity_score <- function(task, learners_trt, learners_cens, control, progress_bar) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_propensity_score(task, fold, learners_trt, learners_cens, control, progress_bar)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  propensity_scores <- lapply(seq_len(task$time_horizon), function(time) {
    this_time <- lapply(seq_along(task$folds), function(fold) {
      ans[[fold]]$propensity_score[, , time]
    })
    recombine(Reduce(rbind, this_time), task$folds)
  })

  ans <- list(
    propensity_score = simplify2array(propensity_scores),
    prob_observed = recombine(rbind_depth(ans, "prob_observed"), task$folds),
    learner_treatment_summary = rbind_depth(ans, "learner_treatment_summary"),
    learner_cens_summary = rbind_depth(ans, "learner_cens_summary")
  )

  # TODO: trimming
  # ans$propensity_score <- trim(ans$propensity_score, control$.trim)
  ans
}

estimate_propensity_score <- function(task, fold, learners_trt, learners_cens, control, progress_bar) {
  time_horizon <- task$time_horizon
  data <- get_folded_data(task$natural, task$folds, fold)
  levels <- task$support(time_horizon)
  number_levels <- length(levels)

  # Local caching
  valid <- data$valid
  train <- data$train
  censvars <- task$vars$C
  trtvars <- task$vars$A
  idvar <- "..i..lmtp_id"

  # Create arrays to store predictions
  propensity_scores <- array(
    NA_real_,
    c(nrow(valid), number_levels, time_horizon),
    dimnames = list(NULL, levels, NULL)
  )

  prob_observed <- matrix(
    ifelse(is.null(task$vars$C), 1, 0),
    nrow = nrow(valid), ncol = time_horizon
  )

  learner_treatment_summary <- NULL
  learner_cens_summary <- NULL

  # TODO: fix issues with point-treatment survival
  for (time in seq_len(time_horizon)) {
    # Get indices of observations that aren't censored and haven't experienced
    # the outcome or the competing risk
    i <- task$observed(train, time - 1) %and% task$is_at_risk(train, time)
    iv <- task$observed(valid, time - 1) %and% task$is_at_risk(valid, time)

    # Current treatment variable
    this_treatment <- current_trt(trtvars, time)
    # Current censoring variable
    this_censoring <- censvars[time]

    # Covariates
    vars <- c(idvar, task$vars$history("A", time))

    # Treatment levels at this time
    levels <- task$support(time_horizon)
    number_levels <- length(levels)

    # One hot encode the treatment
    ohe <- one_hot_encode(train, this_treatment)

    # Loop over K-1 treatment levels as binomial models
    for (l in 2:number_levels) {
      this_level <- colnames(ohe)[l]
      train_ohe <- train

      train_ohe[[this_treatment]] <- ohe[, l]
      fit <- run_ensemble(
        train_ohe[i, c(vars, this_treatment)],
        this_treatment,
        learners_trt, "binomial", idvar, control$.learners_trt_folds
      )

      propensity_scores[iv, this_level, time] <- predict(fit, valid[iv, vars, drop = FALSE])

      # Add fit summary
      learner_treatment_summary <-
        rbind(learner_treatment_summary, summary(fit, time, fold, level = this_level))
    }

    # Add probability for the reference level
    propensity_scores[iv, colnames(ohe)[1], time] <-
      1 - rowSums(propensity_scores[iv, colnames(ohe)[2:number_levels], time, drop = FALSE])

    # Fit model for censoring if there is censoring
    if (!is.null(this_censoring)) {
      fit <- run_ensemble(
        train[i, c(vars, this_treatment, this_censoring)], this_censoring,
        learners_cens, "binomial", "..i..lmtp_id", control$.learners_trt_folds
      )

      # Add fit summary
      learner_cens_summary <- rbind(learner_cens_summary, summary(fit, time, fold))

      # Predict on validation data
      prob_observed[iv, time] <- predict(fit, valid[iv, c(vars, this_treatment), drop = FALSE])
    }

    # Iterate the progress bar
    progress_bar()
  }

  list(propensity_score = propensity_scores,
       prob_observed = prob_observed,
       learner_treatment_summary = learner_treatment_summary,
       learner_cens_summary = learner_cens_summary)
}
