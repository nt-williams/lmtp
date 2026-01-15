cf_propensity_score <- function(task, learners, control, progress_bar) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_propensity_score(task, fold, learners, control, progress_bar)
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
    prob_observed = recombine(rbind_depth(ans, "probability_observed"), task$folds),
    learner_treatment_summary = rbind_depth(ans, "learner_treatment_summary"),
    learner_cens_summary = rbind_depth(ans, "learner_cens_summary")
  )

  # TODO: trimming
  # ans$propensity_score <- trim(ans$propensity_score, control$.trim)
  ans
}

estimate_propensity_score <- function(task, fold, learners, control, pb) {
  data <- get_folded_data(task$natural, task$folds, fold)
  levels <- task$support(task$time_horizon)
  number_levels <- length(levels)

  # Create arrays to store predictions
  propensity_scores <- array(NA_real_, c(nrow(data$valid), number_levels, task$time_horizon), dimnames = list(NULL, levels, NULL))
  probability_observed <- matrix(1, nrow = nrow(data$valid), ncol = task$time_horizon)

  learner_treatment_summary <- NULL
  learner_cens_summary <- learner_treatment_summary

  for (time in seq_len(task$time_horizon)) {
    # Get indices of observations that aren't censored and haven't experienced
    # the outcome or the competing risk
    i <- task$observed(data$train, time - 1) %and% task$is_at_risk(data$train, time)
    iv <- task$observed(data$valid, time - 1) %and% task$is_at_risk(data$valid, time)

    # Current treatment variable
    this_treatment <- current_trt(task$vars$A, time)
    # Current censoring variable
    this_censoring <- task$vars$C[time]

    # Covariates
    vars <- c("..i..lmtp_id", task$vars$history("A", time))

    # Treatment levels at this time
    levels <- task$support(task$time_horizon)
    number_levels <- length(levels)

    # One hot encode the treatment
    ohe <- model.matrix(~ -1 + as.character(data$train[[this_treatment]]))
    colnames(ohe) <- gsub("as.character\\(data\\$train\\[\\[this_treatment\\]\\]\\)", "", colnames(ohe))

    # Loop over K-1 treatment levels as binomial models
    for (l in 2:number_levels) {
      this_level <- colnames(ohe)[l]
      tmp <- data$train
      tmp[[this_treatment]] <- ohe[, l]
      fit <- run_ensemble(
        tmp[i, c(vars, this_treatment)], this_treatment,
        learners, "binomial", "..i..lmtp_id", control$.learners_trt_folds
      )

      propensity_scores[iv, this_level, time] <- predict(fit, data$valid[iv, vars, drop = FALSE])

      # Add fit summary
      learner_treatment_summary <- rbind(learner_treatment_summary, summary(fit, time, fold, level = this_level))
    }

    # Add probability for the reference level
    propensity_scores[iv, colnames(ohe)[1], time] <- 1 - rowSums(propensity_scores[iv, colnames(ohe)[2:number_levels], time, drop = FALSE])

    # Fit model for censoring if there is censoring
    if (!is.null(this_censoring)) {
      fit <- run_ensemble(
        data[i, c(vars, this_treatment, this_censoring)], this_censoring,
        learners, "binomial", "..i..lmtp_id", control$.learners_trt_folds
      )

      # Add fit summary
      learner_cens_summary <- rbind(learner_cens_summary, summary(fit, time, fold, level = NA_character_))
    }

    if (!is.null(this_censoring)) {
      # TODO: predict probability of being observed
    }

    # Iterate the progress bar
    # progress_bar()
  }

  list(propensity_score = propensity_scores,
       probability_observed = probability_observed,
       learner_treatment_summary = learner_treatment_summary,
       learner_cens_summary = learner_cens_summary)
}
