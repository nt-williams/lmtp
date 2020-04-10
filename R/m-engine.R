
# the engine for the initial estimator of m through super learner
estimate_sub <- function(data, shifted, outcome, node_list, C,
                         tau, outcome_type, learner_stack = NULL,
                         m, pb) {

  if (tau > 0) {

    # setup
    i          <- create_censoring_indicators(data, C, tau)$i
    j          <- create_censoring_indicators(data, C, tau)$j
    pseudo     <- paste0("m", tau)
    fit_task   <- initiate_sl3_task(data[i, ], outcome, node_list[[tau]], outcome_type)
    shift_task <- suppressWarnings(initiate_sl3_task(shifted[j, ], outcome, node_list[[tau]], outcome_type))
    ensemble   <- initiate_ensemble(outcome_type, learner_stack)

    # progress bar
    progress_progress_bar(pb)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on shifted data
    m[j, tau]            <-
      shifted[j, pseudo] <-
      data[j, pseudo]    <-
      bound(predict_sl3(fit, shift_task))

    # recursion
    estimate_sub(data = data,
                 shifted = shifted,
                 outcome = pseudo,
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
estimate_tmle <- function(data, shifted, outcome, node_list, C, tau, max,
                          outcome_type, m_natural, m_shifted, r,
                          learner_stack = NULL, pb) {

  if (tau > 0) {

    # setup
    i             <- create_censoring_indicators(data, C, tau)$i
    j             <- create_censoring_indicators(data, C, tau)$j
    pseudo        <- paste0("m", tau)
    fit_task      <- initiate_sl3_task(data[i, ], outcome, node_list[[tau]], outcome_type)
    no_shift_task <- suppressWarnings(initiate_sl3_task(data[j, ], outcome, node_list[[tau]], outcome_type))
    shift_task    <- suppressWarnings(initiate_sl3_task(shifted[j, ], outcome, node_list[[tau]], outcome_type))
    ensemble      <- initiate_ensemble(outcome_type, learner_stack)

    # progress bar
    progress_progress_bar(pb)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on data
    m_natural[j, tau] <- bound(predict_sl3(fit, no_shift_task))
    m_shifted[j, tau] <- bound(predict_sl3(fit, shift_task))

    # tilt estimates
    fit <- suppressWarnings(glm(data[i, outcome] ~ offset(qlogis(m_natural[i, tau])), weights = r[i, tau], family = "binomial"))

    # updating the unshifted estimate
    m_natural[, tau] <- bound(plogis(qlogis(m_natural[, tau]) + coef(fit)))

    # updating the shifted estimate
    m_shifted[, tau]    <-
      shifted[, pseudo] <-
      data[, pseudo]    <-
      bound(plogis(qlogis(m_shifted[, tau]) + coef(fit)))

    # recursion
    estimate_tmle(data = data,
                  shifted = shifted,
                  outcome = pseudo,
                  node_list = node_list,
                  C = C,
                  tau = tau - 1,
                  max = max,
                  outcome_type = "quasibinomial",
                  m_natural = m_natural,
                  m_shifted = m_shifted,
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
estimate_sdr <- function(data, shifted, outcome, node_list, C,
                         tau, max, outcome_type, learner_stack = NULL,
                         m_shifted, m_natural, r, pb) {

  if (tau > 0) {

    # global setup
    pseudo <- paste0("m", tau)

    # progress bar
    progress_progress_bar(pb)

    if (tau == max) {

      # setup
      i             <- create_censoring_indicators(data, C, tau)$i
      j             <- create_censoring_indicators(data, C, tau)$j
      fit_task      <- initiate_sl3_task(data[i, ], outcome, node_list[[tau]], outcome_type)
      no_shift_task <- suppressWarnings(initiate_sl3_task(data[j, ], outcome, node_list[[tau]], outcome_type))
      shift_task    <- suppressWarnings(initiate_sl3_task(shifted[j, ], outcome, node_list[[tau]], outcome_type))
      ensemble      <- initiate_ensemble(outcome_type, learner_stack)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # predict
      m_natural[j, tau]    <- bound(predict_sl3(fit, no_shift_task))
      m_shifted[j, tau]    <-
        shifted[j, pseudo] <-
        data[j, pseudo]    <-
        bound(predict_sl3(fit, shift_task))

    } else if (tau < max) {

      # setup
      i <- create_censoring_indicators(data, C, tau + 1)$i
      j <- create_censoring_indicators(data, C, tau)$j

      # outcome transformation
      z              <- use_dens_ratio(r, tau, NULL, max, "sdr")
      ms             <- m_shifted[, (tau + 2):(max + 1), drop = FALSE] - m_natural[, (tau + 1):max, drop = FALSE]
      mt             <- m_shifted[, tau + 1]
      data[, pseudo] <- shifted[, pseudo] <- rowSums(z * ms) + mt

      # run SL on outcome transformation and get predictions
      fit_task      <- initiate_sl3_task(data[i, ], pseudo, node_list[[tau]], outcome_type)
      no_shift_task <- suppressWarnings(initiate_sl3_task(data[j, ], outcome, node_list[[tau]], outcome_type))
      shift_task    <- suppressWarnings(initiate_sl3_task(shifted[j, ], outcome, node_list[[tau]], outcome_type))
      ensemble      <- initiate_ensemble(outcome_type, learner_stack)
      fit           <- run_ensemble(ensemble, fit_task)

      # predictions
      m_natural[j, tau] <- predict_sl3(fit, no_shift_task)
      m_shifted[j, tau] <- predict_sl3(fit, shift_task)

    }

    # recursion
    estimate_sdr(data = data,
                 shifted = shifted,
                 outcome = pseudo,
                 node_list = node_list,
                 C = C,
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
