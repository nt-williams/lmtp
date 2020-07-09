
estimate_sub <- function(training, shifted, validation, outcome,
                         node_list, C, deterministic, tau, outcome_type,
                         learners = NULL, m, pb, sl_weights) {

  if (tau > 0) {
    # setup
    i          <- create_censoring_indicators(training, C, tau)$i
    js         <- create_censoring_indicators(shifted, C, tau)$j
    jv         <- create_censoring_indicators(validation, C, tau)$j
    dt         <- create_determ_indicators(training, deterministic, tau)
    dv         <- create_determ_indicators(validation, deterministic, tau)
    pseudo     <- paste0("psi", tau)
    fit_task   <- initiate_sl3_task(training[i & !dt, ], outcome, node_list[[tau]], outcome_type, "lmtp_id")
    shift_task <- sw(initiate_sl3_task(shifted[js, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
    valid_task <- sw(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
    ensemble   <- initiate_ensemble(outcome_type, check_variation(training[i, ], outcome, learners))

    # progress bar
    pb()

    # run SL
    fit <- run_ensemble(ensemble, fit_task, envir = environment())
    sl_weights[[tau]] <- extract_sl_weights(fit)

    # predict on shifted data for training
    training[js, pseudo] <- bound(predict_sl3(fit, shift_task, envir = environment()))
    training[dt, pseudo] <- 1

    # predict on validation shifted data
    m[jv, tau] <- bound(predict_sl3(fit, valid_task, envir = environment()))
    m[dv, tau] <- 1

    # recursion
    estimate_sub(training = training,
                 shifted = shifted,
                 validation = validation,
                 outcome = pseudo,
                 node_list = node_list,
                 C = C,
                 deterministic = deterministic,
                 tau = tau - 1,
                 outcome_type = "continuous",
                 learners,
                 m = m,
                 pb = pb,
                 sl_weights = sl_weights)

  } else {
    out <- list(m = m,
                sl_weights = sl_weights)
    # returns
    return(out)
  }
}

estimate_tmle <- function(training, shifted, validation, validation_shifted,
                          outcome, node_list, C, deterministic, tau, outcome_type,
                          m_natural, m_shifted, r, learners = NULL, pb, sl_weights) {

  if (tau > 0) {
    # setup
    i            <- create_censoring_indicators(training, C, tau)$i
    jt           <- create_censoring_indicators(training, C, tau)$j
    jv           <- create_censoring_indicators(validation, C, tau)$j
    dt           <- create_determ_indicators(training, deterministic, tau)
    dv           <- create_determ_indicators(validation, deterministic, tau)
    pseudo       <- paste0("psi", tau)
    fit_task     <- initiate_sl3_task(training[i & !dt, ], outcome, node_list[[tau]], outcome_type, "lmtp_id")
    nshift_task  <- sw(initiate_sl3_task(training[jt, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
    shift_task   <- sw(initiate_sl3_task(shifted[jt, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
    vnshift_task <- sw(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
    vshift_task  <- sw(initiate_sl3_task(validation_shifted[jv, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
    ensemble     <- initiate_ensemble(outcome_type, check_variation(training[i & !dt, ], outcome, learners))

    # progress bar
    pb()

    # run SL
    fit <- run_ensemble(ensemble, fit_task, envir = environment())
    sl_weights[[tau]] <- extract_sl_weights(fit)

    # predict on data
    m_natural$train[jt, tau] <- bound(predict_sl3(fit, nshift_task, envir = environment()))
    m_shifted$train[jt, tau] <- bound(predict_sl3(fit, shift_task, envir = environment()))
    m_natural$valid[jv, tau] <- bound(predict_sl3(fit, vnshift_task, envir = environment()))
    m_shifted$valid[jv, tau] <- bound(predict_sl3(fit, vshift_task, envir = environment()))

    # tilt estimates
    fit <- sw(glm(training[i & !dt, ][[outcome]] ~ offset(qlogis(m_natural$train[i & !dt, tau])),
                                weights = r$train[i & !dt, tau], family = "binomial"))

    # update training estimates
    training[, pseudo] <- bound(plogis(qlogis(m_shifted$train[, tau]) + coef(fit)))
    training[dt, pseudo] <- 1

    # update validation estiamtes
    m_natural$valid[, tau] <- bound(plogis(qlogis(m_natural$valid[, tau]) + coef(fit)))
    m_shifted$valid[, tau] <- bound(plogis(qlogis(m_shifted$valid[, tau]) + coef(fit)))
    m_natural$valid[dv, tau] <- 1
    m_shifted$valid[dv, tau] <- 1

    # recursion
    estimate_tmle(training = training,
                  shifted = shifted,
                  validation = validation,
                  validation_shifted = validation_shifted,
                  outcome = pseudo,
                  node_list = node_list,
                  C = C,
                  deterministic = deterministic,
                  tau = tau - 1,
                  outcome_type = "continuous",
                  m_natural = m_natural,
                  m_shifted = m_shifted,
                  r = r,
                  learners = learners,
                  pb = pb,
                  sl_weights = sl_weights)
  } else {
    # returns
    out <- list(natural = m_natural$valid,
                shifted = m_shifted$valid,
                sl_weights = sl_weights)

    return(out)
  }

}

estimate_sdr <- function(training, shifted, validation, validation_shifted,
                         outcome, node_list, C, deterministic, tau, max, outcome_type,
                         learners = NULL, m_shifted, m_natural, r, pb, sl_weights) {

  if (tau > 0) {

    # global setup
    i      <- create_censoring_indicators(training, C, tau)$i
    jt     <- create_censoring_indicators(training, C, tau)$j
    jv     <- create_censoring_indicators(validation, C, tau)$j
    dt     <- create_determ_indicators(training, deterministic, tau)
    dv     <- create_determ_indicators(validation, deterministic, tau)
    pseudo <- paste0("psi", tau + 1)

    # progress bar
    pb()

    if (tau == max) {

      # setup
      fit_task     <- initiate_sl3_task(training[i & !dt, ], outcome, node_list[[tau]], outcome_type, "lmtp_id")
      nshift_task  <- sw(initiate_sl3_task(training[jt, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
      shift_task   <- sw(initiate_sl3_task(shifted[jt, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
      vnshift_task <- sw(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
      vshift_task  <- sw(initiate_sl3_task(validation_shifted[jv, ], NULL, node_list[[tau]], NULL, "lmtp_id"))
      ensemble     <- initiate_ensemble(outcome_type, check_variation(training[i & !dt, ], outcome, learners))

      # run SL
      fit <- run_ensemble(ensemble, fit_task, envir = environment())
      sl_weights[[tau]] <- extract_sl_weights(fit)

      # predict on training data
      m_natural$train[jt, tau] <- bound(predict_sl3(fit, nshift_task, envir = environment()))
      m_shifted$train[jt, tau] <- bound(predict_sl3(fit, shift_task, envir = environment()))
      m_natural$train[dt, tau] <- 1
      m_shifted$train[dt, tau] <- 1

      # predict on validation data
      m_natural$valid[jv, tau] <- bound(predict_sl3(fit, vnshift_task, envir = environment()))
      m_shifted$valid[jv, tau] <- bound(predict_sl3(fit, vshift_task, envir = environment()))
      m_natural$valid[dv, tau] <- 1
      m_shifted$valid[dv, tau] <- 1

    } else if (tau < max) {

      # outcome transformation
      training[, pseudo]  <-
        shifted[, pseudo] <-
        transform_sdr(ratio_sdr(r$train, tau, max),
                      tau, max, m_shifted$train, m_natural$train)

      # run SL on outcome transformation and get predictions
      fit_task     <- initiate_sl3_task(training[i & !dt, ], pseudo, node_list[[tau]], outcome_type, "lmtp_id")
      nshift_task  <- sw(initiate_sl3_task(training[jt, ], NULL, node_list[[tau]], outcome_type, "lmtp_id"))
      shift_task   <- sw(initiate_sl3_task(shifted[jt, ], NULL, node_list[[tau]], outcome_type, "lmtp_id"))
      vnshift_task <- sw(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], outcome_type, "lmtp_id"))
      vshift_task  <- sw(initiate_sl3_task(validation_shifted[jv, ], NULL, node_list[[tau]], outcome_type, "lmtp_id"))
      ensemble     <- initiate_ensemble(outcome_type, check_variation(training[i & !dt, ], pseudo, learners))

      # run SL
      fit <- run_ensemble(ensemble, fit_task, envir = environment())
      sl_weights[[tau]] <- extract_sl_weights(fit)

      # predictions
      m_natural$train[jt, tau] <- bound(predict_sl3(fit, nshift_task, envir = environment()))
      m_shifted$train[jt, tau] <- bound(predict_sl3(fit, shift_task, envir = environment()))
      m_natural$train[dt, tau] <- 1
      m_shifted$train[dt, tau] <- 1
      m_natural$valid[jv, tau] <- bound(predict_sl3(fit, vnshift_task, envir = environment()))
      m_shifted$valid[jv, tau] <- bound(predict_sl3(fit, vshift_task, envir = environment()))
      m_natural$valid[dv, tau] <- 1
      m_shifted$valid[dv, tau] <- 1

    }

    # recursion
    estimate_sdr(training = training,
                 shifted = shifted,
                 validation = validation,
                 validation_shifted = validation_shifted,
                 outcome = pseudo,
                 node_list = node_list,
                 C = C,
                 deterministic = deterministic,
                 tau = tau - 1,
                 max = max,
                 outcome_type = "continuous",
                 learners = learners,
                 m_shifted = m_shifted,
                 m_natural = m_natural,
                 r = r,
                 pb = pb,
                 sl_weights = sl_weights)
  } else {
    # returns
    out <- list(natural = m_natural$valid,
                shifted = m_shifted$valid,
                sl_weights = sl_weights)

    return(out)
  }

}
