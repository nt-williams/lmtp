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


cf_sequential_regression <- function(task, final_time, scores, learners, bounds = NULL, control, pb = NULL) {
  
  ans <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_flip_sdr(task, final_time, fold, scores, learners, bounds, control, pb)
    },
    seed = TRUE)
  }
  
  ans <- future::value(ans)

  return(
    list(
      seq_reg_ests = abind::abind(lapply(ans, `[[`, "seq_reg_ests"), along = 1)[ 
        order(unlist(lapply(task$folds, `[[`, "validation_set"))), ,, drop = F],
      fits = lapply(ans, function(x) x[["fits"]])
    )
  )
}

estimate_flip_sdr <- function(task, final_time, fold, scores, learners, bounds, control, pb = NULL) {

  # Assign data to training / validation
  data <- get_folded_data(task$natural, task$folds, fold)

  # Get q_scores; phi_scores; ratios into train & validation?
  q_scores <- get_folded_data(scores$q_scores, task$folds, fold)$train
  ratios <- get_folded_data(scores$ratios, task$folds, fold)$train
  phi_scores <- scores$phi_scores[task$folds[[fold]]$training_set, , , drop = FALSE]
  
  # Initialize array for sequential regression estimates in-training
  m_train <- array(dim = c(nrow(data$train), final_time, 2)) 

  # Use m_valid to keep final output sequential regression estimates
  m_valid <- array(dim = c(nrow = nrow(data$valid), final_time, ncol = 2))

  # --- Run sequential regression on training data --- #
  fits <- vector("list", length = final_time)
  for (time in final_time:1) {
    
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
                        outcome_type = ifelse(time != final_time, "continuous", task$outcome_type),
                        id = "..i..lmtp_id",
                        folds = control$.learners_outcome_folds)
    
    # --- Predict under treatment and control ---- # 
    
    # Predict under A_t = 0 and A_t = 1 for all relevant subjects
    # In training data
    untreated_train <- treated_train <- data$train; 
    untreated_train[, task$vars$A[[time]]] <- 0
    treated_train[, task$vars$A[[time]]] <- 1
    
    m_train[i, time, 1] <- predict(fit, untreated_train[i, ], NULL)
    m_train[i, time, 2] <- predict(fit, treated_train[i, ], NULL)
    
    # Over-write for censored or already dead
    m_train[which(!y1), time, ] <- 0 
    m_train[which(!d0), time, ] <- 1
    
    # In validation data
    # More censoring logic: check whether censored in previous timepoint 
    cp1 <- task$observed(data$train, time-1) # censoring in the past = 1
    y1v <- task$at_risk_N(data$valid, time-1)
    d0v <- task$at_risk_D(data$valid, time-1)
    cp1v <- task$observed(data$valid, time-1)
    iv <- ii(cp1v, y1v & d0v)
    
    untreated <- treated <- data$valid; 
    untreated[, task$vars$A[[time]]] <- 0; treated[, task$vars$A[[time]]] <- 1
    
    m_valid[iv, time, 1] <- predict(fit, untreated[iv, ], NULL)
    m_valid[iv, time, 2] <- predict(fit, treated[iv, ], NULL)
    
    # Over-write for censored or already dead
    m_valid[which(!y1v), time, ] <- 0
    m_valid[which(!d0v), time, ] <- 1
    
    # Return fits and estimates
    fits[[time]] <- if (control$.return_full_fits) fit else extract_sl_weights(fit)

    # Construct pseudo outcomes for next round
    data$train$pseudo <- eif_Q(task = task, 
                               data = data$train,
                               seq_ests = m_train,
                               q_scores = q_scores,
                               phi_scores = phi_scores,
                               ratios = ratios,
                               time = time,
                               final_time = final_time)
    
    if(!is.null(bounds)) { # Constrain pseudo-outcome to bounded range if applicable
      data$train$pseudo <- pmax(bounds[[1]], pmin(bounds[[2]], data$train$pseudo))
    }
  }
  
  if (!is.null(pb)) pb()
  
  # Return sequential regression estimates from 
  return(list(seq_reg_ests = m_valid, fits = fits))
}
