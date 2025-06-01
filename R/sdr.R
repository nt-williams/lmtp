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
       fits = lapply(ans, function(x) x[["fits"]]))
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
    y1 <- task$at_risk_N(natural$train, t-1)
    d0 <- task$at_risk_D(natural$train, t-1)
    c1 <- task$observed(natural$train, t)
    i <- ii(c1, y1 & d0)

    history <- task$vars$history("L", t + 1)
    vars <- c("..i..lmtp_id", history, task$vars$Y)

    fit <- run_ensemble(natural$train[i, vars], task$vars$Y,
                        learners,
                        ifelse(t != task$tau, "continuous", task$outcome_type),
                        "..i..lmtp_id",
                        control$.learners_outcome_folds)

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

    cp1 <- task$observed(natural$train, t-1) # censoring in the past = 1
    y1v <- task$at_risk_N(natural$valid, t-1)
    d0v <- task$at_risk_D(natural$valid, t-1)
    cp1v <- task$observed(natural$valid, t-1)

    i <- ii(cp1, y1 & d0)
    iv <- ii(cp1v, y1v & d0v)

    under_shift_train <- natural$train[i, c("..i..lmtp_id", history)]
    under_shift_train[, A_t] <- shifted$train[i, A_t]

    m_natural_train[i, t] <- predict(fit, natural$train[i, ], NULL)
    m_shifted_train[i, t] <- predict(fit, under_shift_train, NULL)

    m_natural_train[which(!y1), t] <- 0
    m_natural_train[which(!d0), t] <- 1
    m_shifted_train[which(!y1), t] <- 0
    m_shifted_train[which(!d0), t] <- 1

    under_shift_valid <- natural$valid[iv, c("..i..lmtp_id", history)]
    under_shift_valid[, A_t] <- shifted$valid[iv, A_t]

    m_natural_valid[iv, t] <- predict(fit, natural$valid[iv, ], NULL)
    m_shifted_valid[iv, t] <- predict(fit, under_shift_valid, NULL)

    m_natural_valid[which(!y1v), t] <- 0
    m_natural_valid[which(!d0v), t] <- 1
    m_shifted_valid[which(!y1v), t] <- 0
    m_shifted_valid[which(!d0v), t] <- 1

    natural$train[, task$vars$Y] <- eif(ratios, m_shifted_train, m_natural_train, t)

    pb()
  }

  list(natural = m_natural_valid,
       shifted = m_shifted_valid,
       fits = fits)
}




cf_sequential_regression <- function(task, learners, control, time, pb = NULL) {
  
  folds <- task$folds
  nf <- length(folds)
  ans <- lapply(seq_len(nf), function(f) estimate_flip_sdr(task, f, learners, control, time, pb))
  
  return(list(seq_reg_ests = recombine(rbind_depth(ans, "seq_reg_ests"), task$folds),
              fits = lapply(ans, function(x) x[["fits"]])))
}

estimate_flip_sdr <- function(task, fold, learners, control, time, pb = NULL) {
  
  # --- Run sequential regression on training data --- #
  
  # Assign data to training / validation
  data <- get_folded_data(task$natural, task$folds, fold)
  m_valid <- matrix(nrow = nrow(data$valid), ncol = 2)
  
  # Logic for dealing with censoring and competing events
  y1 <- task$at_risk_N(data$train, time-1)
  d0 <- task$at_risk_D(data$train, time-1)
  c1 <- task$observed(data$train, time)
  i <- ii(c1, y1 & d0) # i = index of subjects to use for regression
  
  # Get predictors for sequential regression
  history <- task$vars$history("L", time + 1)
  vars <- c("..i..lmtp_id", history, "pseudo")
  
  # Run sequential regression to estimate m_t(A_t, H_t)
  fit <- run_ensemble(data = data$train[i, vars], 
                      y = "pseudo", # Repeatedly updated pseudo-outcome
                      learners = learners,
                      outcome_type = ifelse(time != task$tau, "continuous", task$outcome_type),
                      id = "..i..lmtp_id",
                      folds = control$.learners_outcome_folds)
  
  # --- Predict under treatment and control ---- # 
  
  # More censoring logic: check whether censored in previous timepoint 
  cp1 <- task$observed(data$train, time-1) # censoring in the past = 1
  y1v <- task$at_risk_N(data$valid, time-1)
  d0v <- task$at_risk_D(data$valid, time-1)
  cp1v <- task$observed(data$valid, time-1)
  iv <- ii(cp1v, y1v & d0v)
  
  # Predict under A_t = 0 and A_t = 1 for all relevant subjects
  untreated <- treated <- data$valid; 
  untreated[, task$vars$A[[time]]] <- 0; treated[, task$vars$A[[time]]] <- 1
  
  m_valid[iv, 1] <- predict(fit, untreated[iv, ], NULL)
  m_valid[iv, 2] <- predict(fit, treated[iv, ], NULL)
  
  # Over-write for censored or already dead
  m_valid[which(!y1v), ] <- 0
  m_valid[which(!d0v), ] <- 1
  
  if (!is.null(pb)) pb()
  return(list(seq_reg_ests = m_valid, fits = if (control$.return_full_fits) fit else extract_sl_weights(fit)))
}
