
# the engine for the initial estimator of m through super learner
estimate_sub <- function(data, shifted, Y, node_list, C,
                         tau, outcome_type, learner_stack = NULL,
                         m, pb) {

  if (tau > 0) {

    # setup
    data <- data[data[[C[tau]]] == 1, ]
    shifted <- data[data[[C[tau]]] == 1, ]
    fit_task <- initiate_sl3_task(data, Y, node_list[[tau]], outcome_type)
    pred_task <- initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type)
    ensemble <- initiate_ensemble(outcome_type, learner_stack)

    # progress bar
    progress_progress_bar(pb)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on shifted data
    pseudo <- paste0("m", tau)
    m[, tau] <- shifted[, pseudo] <- data[, pseudo] <- bound(predict_sl3(fit, pred_task))

    # recursion
    estimate_sub(data = data,
                 shifted = shifted,
                 Y = pseudo,
                 node_list = node_list,
                 C = C,
                 tau = tau - 1,
                 outcome_type = "quasibinomial",
                 learner_stack,
                 m = m,
                 pb = pb)

  } else {
    # returns
    return(m)
  }
}

# the engine for the TML estimator
estimate_tmle <- function(data, shifted, Y, node_list, C, tau, max,
                          outcome_type, m_natural, m_shifted,
                          m_natural_initial, m_shifted_initial, r,
                          learner_stack = NULL, pb) {

  if (tau > 0) {

    # setup
    i <- data[[C[tau]]] == 1
    fit_task <- initiate_sl3_task(data[i, ], Y, node_list[[tau]], outcome_type)
    no_shift_task <- suppressWarnings(initiate_sl3_task(data, Y, node_list[[tau]], outcome_type))
    shift_task <- suppressWarnings(initiate_sl3_task(shifted, Y, node_list[[tau]], outcome_type))
    ensemble <- initiate_ensemble(outcome_type, learner_stack)

    # progress bar
    progress_progress_bar(pb)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on data
    pseudo <- paste0("m", tau)
    m_natural_initial[, tau] <- bound(predict_sl3(fit, no_shift_task))
    m_shifted_initial[, tau] <- bound(predict_sl3(fit, shift_task))

    # tilt estimates
    fit <- suppressWarnings(glm(data[i, Y] ~ offset(qlogis(m_natural_initial[i, tau])), weights = r[i, tau], family = "binomial"))

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
                  C = C,
                  tau = tau - 1,
                  max = max,
                  outcome_type = "quasibinomial",
                  m_natural = m_natural,
                  m_shifted = m_shifted,
                  m_natural_initial = m_natural_initial,
                  m_shifted_initial = m_shifted_initial,
                  r = r,
                  learner_stack = learner_stack,
                  pb = pb)
  } else {
    # returns
    out <- list(natural = m_natural,
                shifted = m_shifted)

    return(out)
  }

}

# the engine for the sdr estimator
estimate_sdr <- function(data, shifted, Y, node_list,
                         tau, max, outcome_type, learner_stack = NULL,
                         m_shifted, m_natural, r, pb) {

  if (tau > 0) {

    # progress bar
    progress_progress_bar(pb)

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
      z <- use_dens_ratio(r, tau, NULL, max, "sdr")
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
                 r = r,
                 pb = pb)
  } else {
    # returns
    out <- list(natural = m_natural,
                shifted = m_shifted)

    return(out)
  }

}
