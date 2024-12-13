cf_sdr <- function(task, ratios, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_sdr(task, fold, ratios, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
       shifted = recombine(rbind_depth(ans, "shifted"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_sdr <- function(task, fold, ratios, learners, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  ratios <- get_folded_data(ratios, task$folds, fold)$train

  m_natural_train <- matrix(nrow = nrow(natural$train), ncol = task$tau + 1)
  m_shifted_train <- m_natural_train

  m_natural_valid <- matrix(nrow = nrow(natural$valid), ncol = task$tau + 1)
  m_shifted_valid <- m_natural_valid

  m_shifted_train[, task$tau + 1] <- natural$train[[task$vars$Y]]
  m_shifted_valid[, task$tau + 1] <- natural$valid[[task$vars$Y]]

  fits <- vector("list", length = task$tau)
  for (t in task$tau:1) {
    i <- ii(task$observed(natural$train, t), task$at_risk(natural$train, t))

    history <- task$vars$history("L", t + 1)
    vars <- c("..i..lmtp_id", history, task$vars$Y)

    fit <- run_ensemble(natural$train[i, vars], task$vars$Y,
                        learners,
                        ifelse(t != task$tau, "continuous", task$outcome_type),
                        "..i..lmtp_id",
                        control$.learners_outcome_folds,
                        control$.discrete,
                        control$.info)

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    if (length(task$vars$A) > 1) {
      A_t <- task$vars$A[[t]]
    } else {
      A_t <- task$vars$A[[1]]
    }

    i <- ii(task$observed(natural$train, t - 1), task$at_risk(natural$train, t))
    under_shift_train <- natural$train[i, c("..i..lmtp_id", history)]
    under_shift_train[, A_t] <- shifted$train[i, A_t]

    m_natural_train[i, t] <- predict(fit, natural$train[i, ], NULL)
    m_shifted_train[i, t] <- predict(fit, under_shift_train, NULL)

    i <- ii(task$observed(natural$valid, t - 1), task$at_risk(natural$valid, t))
    under_shift_valid <- natural$valid[i, c("..i..lmtp_id", history)]
    under_shift_valid[, A_t] <- shifted$valid[i, A_t]

    m_natural_valid[i, t] <- predict(fit, natural$valid[i, ], NULL)
    m_shifted_valid[i, t] <- predict(fit, under_shift_valid, NULL)

    m_natural_valid[which(!(task$at_risk(natural$valid, t))), t] <- 0
    m_shifted_valid[which(!(task$at_risk(natural$valid, t))), t] <- 0
    m_natural_train[which(!(task$at_risk(natural$train, t))), t] <- 0
    m_shifted_train[which(!(task$at_risk(natural$train, t))), t] <- 0

    natural$train[, task$vars$Y] <- eif(ratios, m_shifted_train, m_natural_train, t)
    natural$train[which(!(task$at_risk(natural$train, t))), task$vars$Y] <- 0

    pb()
  }

  list(natural = m_natural_valid,
       shifted = m_shifted_valid,
       fits = fits)
}
