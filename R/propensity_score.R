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

  # # TODO: trimming
  # browser()
  # ans$propensity_score <- trim(ans$propensity_score, control$.trim)
  ans
}

estimate_propensity_score <- function(task, fold, learners_trt, learners_cens, control, pb) {
  time_horizon <- task$time_horizon
  data <- get_folded_data(task$natural, task$folds, fold)
  levels <- task$support(1)
  number_levels <- length(levels)

  # Local caching
  valid <- data$valid
  train <- data$train
  censvars <- task$vars$C
  trtvars <- task$vars$A
  idvar <- "..i..lmtp_id"

  # Create arrays to store predictions
  propensity_scores <- array(NA_real_, c(nrow(valid), number_levels, time_horizon),
                             dimnames = list(NULL, levels, NULL))
  prob_observed <- matrix(ifelse(is.null(task$vars$C), 1, NA_real_),
                          nrow = nrow(valid), ncol = time_horizon)

  learner_treatment_summary <- NULL
  learner_cens_summary <- NULL

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
    levels <- task$support(time)
    number_levels <- length(levels)
    
    # Loop over K-1 treatment levels as binomial models
    # Need to handle degenerate case of no variation in treatment
    if (number_levels > 1) {
    # One hot encode the treatment
    ohe <- one_hot_encode(train[i, ], this_treatment)
      for (l in 2:number_levels) {
        this_level <- colnames(ohe)[l]
        fit <- run_ensemble(
          collapse::ftransform(train[i, ], this_treatment = ohe[, l])[, c(vars, this_treatment)],
          this_treatment,
          learners_trt, "binomial", idvar, control$.learners_trt_folds,
          control$.discrete, control$.info
        )

        propensity_scores[iv, this_level, time] <- predict(fit, valid[iv, vars, drop = FALSE])
      }

      # Add probability for the reference level
      propensity_scores[iv, colnames(ohe)[1], time] <-
        1 - rowSums(propensity_scores[iv, colnames(ohe)[2:number_levels], time, drop = FALSE])
    } else {
      this_level <- as.character(unique(train[[this_treatment]]))
      fit <- run_ensemble(
        collapse::ftransform(train[i, ], this_treatment = train[[this_treatment]])[, c(vars, this_treatment)],
        this_treatment,
        "mean", "binomial", idvar, control$.learners_trt_folds,
        control$.discrete, control$.info
      )

      propensity_scores[iv, this_level, time] <- predict(fit, valid[iv, vars, drop = FALSE])

      # Add probability for the reference level
      propensity_scores[iv, colnames(ohe)[1], time] <-
        1 - rowSums(propensity_scores[iv, this_level, time, drop = FALSE])
    }

    # Add fit summary
    learner_treatment_summary <-
      rbind(learner_treatment_summary, summary(fit, time, fold, level = this_level))

    # Fit model for censoring if there is censoring
    if (!is.null(this_censoring)) {
      fit <- run_ensemble(
        train[i, c(task$vars$history_cens(time), this_censoring)], this_censoring,
        learners_cens, "binomial", "..i..lmtp_id", control$.learners_trt_folds,
        control$.discrete, control$.info
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
