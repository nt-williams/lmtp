
# the engine for the parametric substitution estimator
estimate_m_glm <- function(data, shifted, Y,
                           node_list, tau, family, m) {

  if (tau > 0) {
    # setup
    form <- as.formula(paste(Y, paste(node_list[[tau]], collapse = " + "), sep = " ~ "))

    # fit glm
    fit <- glm(form, data = data, family = family)

    # predict on shifted data
    pseudo <- paste0("m", tau)
    data[, pseudo] <- m[, tau] <- qlogis(bound(plogis(predict(fit, newdata = shifted))))

    # recursion
    estimate_m_glm(data = data,
                   shifted = shifted,
                   tau = tau - 1,
                   node_list = node_list,
                   Y = pseudo,
                   m = m,
                   family = "gaussian")
  } else {
    # when t = 1 return matrix m
    return(m)
  }
}

# the engine for the sdr estimator
estimate_sdr <- function(data, shifted, Y, node_list,
                         tau, max, outcome_type, learner_stack = NULL,
                         m_shifted, m_natural, r) {

  if (tau > 0) {

    if (tau == max) {

      # setup
      pseudo <- paste0("m", tau)
      fit_task <- initiate_sl3_task(data, Y, node_list[[tau]], outcome_type)
      pred_task <- initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type)
      ensemble <- initiate_ensemble(outcome_type, learner_stack)
      m_natural <- cbind(m_natural, data[, Y])
      m_shifted <- cbind(m_shifted, data[, Y])

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # predict on shifted data
      m_natural[, tau] <- bound(predict_sl3(fit, fit_task))
      m_shifted[, tau] <- shifted[, pseudo] <- data[, pseudo] <- bound(predict_sl3(fit, pred_task))

    } else if (tau < max) {

      # setup
      pseudo <- "y_sdr"

      # outcome transformation
      z <- use_r_sdr(r, tau, max)
      ms <- m_shifted[, (tau + 2):(max + 1), drop = FALSE] - m_natural[, (tau + 1):max, drop = FALSE]
      mt <- m_shifted[, tau + 1]
      data[, pseudo] <- shifted[, pseudo] <- rowSums(z * ms) + mt

      # run SL on outcome transformation and get predictions
      fit_task <- initiate_sl3_task(data, pseudo, node_list[[tau]], outcome_type)
      pred_task <- initiate_sl3_task(shifted, pseudo, node_list[[tau]], outcome_type)
      ensemble <- initiate_ensemble(outcome_type, learner_stack)
      fit <- run_ensemble(ensemble, fit_task)
      m_natural[, tau] <- predict_sl3(fit, fit_task)
      m_shifted[, tau] <- predict_sl3(fit, pred_task)

    }

    # recursion
    estimate_sdr(data = data,
                 shifted = shifted,
                 Y = pseudo,
                 node_list = node_list,
                 tau = tau - 1,
                 max = max,
                 outcome_type = "continuous",
                 learner_stack = learner_stack,
                 m_shifted = m_shifted,
                 m_natural = m_natural,
                 r = r)
  } else {
    # returns
    return(m_shifted)
  }

}

# the engine for the initial estimator of m through super learner
estimate_m_sl <- function(data, shifted, Y, node_list,
                          tau, outcome_type, learner_stack = NULL,
                          estimator, m) {

  if (tau > 0) {
    # setup
    fit_task <- initiate_sl3_task(data, Y, node_list[[tau]], outcome_type)
    pred_task <- initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type)
    ensemble <- initiate_ensemble(outcome_type, learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on shifted data
    pseudo <- paste0("m", tau)
    m[, tau] <- shifted[, pseudo] <- data[, pseudo] <- bound(predict_sl3(fit, pred_task))

    # recursion
    estimate_m_sl(data = data,
                  shifted = shifted,
                  Y = pseudo,
                  node_list = node_list,
                  tau = tau - 1,
                  outcome_type = "quasibinomial",
                  learner_stack,
                  m = m)

  } else {
    # when t = 1 return matrix m
    return(m)
  }
}

# the engine for the TML estimator
estimate_tmle <- function(data, shifted, Y, node_list, tau,
                          outcome_type, m_natural, m_shifted,
                          m_natural_initial, m_shifted_initial, r,
                          learner_stack = NULL) {

  if (tau > 0) {
    # setup
    fit_task <- initiate_sl3_task(data, Y, node_list[[tau]], outcome_type)
    pred_task <- initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type)
    ensemble <- initiate_ensemble(outcome_type, learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on data
    pseudo <- paste0("m", tau)
    m_natural_initial[, tau] <- bound(predict_sl3(fit, fit_task))
    m_shifted_initial[, tau] <- bound(predict_sl3(fit, pred_task))

    # tilt estimates
    fit <- suppressWarnings(glm(data[, Y] ~ offset(qlogis(m_natural_initial[, tau])),
                                weights = r$rn[, tau], family = "binomial"))

    # updating the unshifted estimate
    m_natural[, tau] <- bound(plogis(qlogis(m_natural_initial[, tau]) + coef(fit)))

    # updating the shifted estimate
    m_shifted[, tau] <-
      shifted[, pseudo] <-
      data[, pseudo] <-
      bound(plogis(qlogis(m_shifted_initial[, tau]) + coef(fit)))

    # recursion
    estimate_tmle(data = data,
                  shifted = shifted,
                  Y = pseudo,
                  node_list = node_list,
                  tau = tau - 1,
                  outcome_type = "quasibinomial",
                  m_natural = m_natural,
                  m_shifted = m_shifted,
                  m_natural_initial = m_natural_initial,
                  m_shifted_initial = m_shifted_initial,
                  r = r,
                  learner_stack = learner_stack)
  } else {
    # returns tilted m when finished
    return(m_shifted)
  }

}
