
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

# the engine for the initial estimator of m through super learner
estimate_m_sl <- function(data, shifted, Y, node_list,
                          tau, outcome_type, learner_stack = NULL, m) {

  if (tau > 0) {
    # setup
    fit_task <- initiate_sl3_task(data, Y, node_list[[tau]], outcome_type)
    pred_task <- initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type)
    ensemble <- initiate_ensemble(outcome_type, learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on shifted data
    pseudo <- paste0("m", tau)
    m[, tau] <- shifted[, pseudo] <- data[, pseudo] <- bound(predict_sl3_nondensity(fit, pred_task))

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
    m_natural_initial[, tau] <- bound(predict_sl3_nondensity(fit, fit_task))
    m_shifted_initial[, tau] <- bound(predict_sl3_nondensity(fit, pred_task))

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
                  outcome_type = "continuous",
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
