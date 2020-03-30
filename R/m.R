
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
                          tau, outcome_type, learners = NULL, m) {

  if (tau > 0) {
    # setup
    ensemble <- initiate_ensemble(data, Y, node_list[[tau]], outcome_type, learners)
    to_predict <- initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type)

    # run SL
    fit <- run_ensemble(ensemble)

    # predict on shifted data
    pseudo <- paste0("m", tau)
    m[, tau] <- shifted[, pseudo] <- data[, pseudo] <- bound(predict_sl3_nondensity(fit, to_predict))

    # recursion
    estimate_m_sl(data = data,
                  shifted = shifted,
                  Y = pseudo,
                  node_list = node_list,
                  tau = tau - 1,
                  outcome_type = "quasibinomial",
                  learners,
                  m = m)

  } else {
    # when t = 1 return matrix m
    return(m)
  }
}
