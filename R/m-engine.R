
# the engine for the initial estimator of m through super learner
#' @export
estimate_sub <- function(training, shifted, validation, outcome, node_list, C,
                         tau, outcome_type, learners = NULL, m, pb) {

  if (tau > 0) {

    # setup
    i          <- create_censoring_indicators(training, C, tau)$i
    js         <- create_censoring_indicators(shifted, C, tau)$j
    jv         <- create_censoring_indicators(validation, C, tau)$j
    pseudo     <- paste0("m", tau)
    fit_task   <- initiate_sl3_task(training[i, ], outcome, node_list[[tau]], outcome_type)
    shift_task <- suppressWarnings(initiate_sl3_task(shifted[js, ], NULL, node_list[[tau]], NULL))
    valid_task <- suppressWarnings(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL))
    ensemble   <- initiate_ensemble(outcome_type, check_variation(training[i, ], outcome, learners))

    # progress bar
    progress_progress_bar(pb)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on shifted data for training
    training[js, pseudo] <- bound(predict_sl3(fit, shift_task))

    # predict on validation shifted data
    m[jv, tau] <- bound(predict_sl3(fit, valid_task))

    # recursion
    estimate_sub(training = training,
                 shifted = shifted,
                 validation = validation,
                 outcome = pseudo,
                 node_list = node_list,
                 C = C,
                 tau = tau - 1,
                 outcome_type = "continuous",
                 learners,
                 m = m,
                 pb = pb)

  } else {
    # returns
    return(m)
  }
}

# the engine for the TML estimator
#' @export
estimate_tmle <- function(training, shifted, validation, validation_shifted,
                          outcome, node_list, C, tau, outcome_type,
                          m_natural, m_shifted, r, learners = NULL, pb) {

  if (tau > 0) {

    # setup
    i            <- create_censoring_indicators(training, C, tau)$i
    jt           <- create_censoring_indicators(training, C, tau)$j
    jv           <- create_censoring_indicators(validation, C, tau)$j
    pseudo       <- paste0("m", tau)
    fit_task     <- initiate_sl3_task(training[i, ], outcome, node_list[[tau]], outcome_type)
    nshift_task  <- suppressWarnings(initiate_sl3_task(training[jt, ], NULL, node_list[[tau]], NULL))
    shift_task   <- suppressWarnings(initiate_sl3_task(shifted[jt, ], NULL, node_list[[tau]], NULL))
    vnshift_task <- suppressWarnings(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL))
    vshift_task  <- suppressWarnings(initiate_sl3_task(validation_shifted[jv, ], NULL, node_list[[tau]], NULL))
    ensemble     <- initiate_ensemble(outcome_type, check_variation(training[i, ], outcome, learners))

    # progress bar
    progress_progress_bar(pb)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # predict on data
    m_natural$train[jt, tau] <- bound(predict_sl3(fit, nshift_task))
    m_shifted$train[jt, tau] <- bound(predict_sl3(fit, shift_task))
    m_natural$valid[jv, tau] <- bound(predict_sl3(fit, vnshift_task))
    m_shifted$valid[jv, tau] <- bound(predict_sl3(fit, vshift_task))

    # tilt estimates
    fit <- suppressWarnings(glm(training[i, outcome] ~ offset(qlogis(m_natural$train[i, tau])),
                                weights = r$train[i, tau], family = "binomial"))

    # update training estimates
    training[, pseudo] <- bound(plogis(qlogis(m_shifted$train[, tau]) + coef(fit)))

    # update validation estiamtes
    m_natural$valid[, tau] <- bound(plogis(qlogis(m_natural$valid[, tau]) + coef(fit)))
    m_shifted$valid[, tau] <- bound(plogis(qlogis(m_shifted$valid[, tau]) + coef(fit)))

    # recursion
    estimate_tmle(training = training,
                  shifted = shifted,
                  validation = validation,
                  validation_shifted = validation_shifted,
                  outcome = pseudo,
                  node_list = node_list,
                  C = C,
                  tau = tau - 1,
                  outcome_type = "continuous",
                  m_natural = m_natural,
                  m_shifted = m_shifted,
                  r = r,
                  learners = learners,
                  pb = pb)
  } else {
    # returns
    out <- list(natural = m_natural$valid,
                shifted = m_shifted$valid)

    return(out)
  }

}

# the engine for the sdr estimator
#' @export
estimate_sdr <- function(training, shifted, validation, validation_shifted,
                         outcome, node_list, C, tau, max, outcome_type,
                         learners = NULL, m_shifted, m_natural, r, pb) {

  if (tau > 0) {

    # global setup
    pseudo <- paste0("m", tau)

    # progress bar
    progress_progress_bar(pb)

    if (tau == max) {

      # setup
      i            <- create_censoring_indicators(training, C, tau)$i
      jt           <- create_censoring_indicators(training, C, tau)$j
      jv           <- create_censoring_indicators(validation, C, tau)$j
      fit_task     <- initiate_sl3_task(training[i, ], outcome, node_list[[tau]], outcome_type)
      nshift_task  <- suppressWarnings(initiate_sl3_task(training[jt, ], NULL, node_list[[tau]], NULL))
      shift_task   <- suppressWarnings(initiate_sl3_task(shifted[jt, ], NULL, node_list[[tau]], NULL))
      vnshift_task <- suppressWarnings(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL))
      vshift_task  <- suppressWarnings(initiate_sl3_task(validation_shifted[jv, ], NULL, node_list[[tau]], NULL))
      ensemble     <- initiate_ensemble(outcome_type, check_variation(training[i, ], outcome, learners))

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # predict on training data
      m_natural$train[jt, tau] <- bound(predict_sl3(fit, nshift_task))

      training[jt, pseudo]       <-
        m_shifted$train[jt, tau] <-
        bound(predict_sl3(fit, shift_task))

      # predict on validation data
      m_natural$valid[jv, tau] <- bound(predict_sl3(fit, vnshift_task))
      m_shifted$valid[jv, tau] <- bound(predict_sl3(fit, vshift_task))

    } else if (tau < max) {

      # setup
      i  <- create_censoring_indicators(training, C, tau + 1)$i
      jt <- create_censoring_indicators(training, C, tau)$j
      jv <- create_censoring_indicators(validation, C, tau)$j

      # outcome transformation
      training[, pseudo]  <-
        shifted[, pseudo] <-
        transform_sdr(ratio_sdr(r$train, tau, max),
                      tau, max, m_shifted$train, m_natural$train)

      # run SL on outcome transformation and get predictions
      fit_task     <- initiate_sl3_task(training[i, ], pseudo, node_list[[tau]], outcome_type)
      nshift_task  <- suppressWarnings(initiate_sl3_task(training[jt, ], NULL, node_list[[tau]], NULL))
      shift_task   <- suppressWarnings(initiate_sl3_task(shifted[jt, ], NULL, node_list[[tau]], NULL))
      vnshift_task <- suppressWarnings(initiate_sl3_task(validation[jv, ], NULL, node_list[[tau]], NULL))
      vshift_task  <- suppressWarnings(initiate_sl3_task(validation_shifted[jv, ], NULL, node_list[[tau]], NULL))
      ensemble     <- initiate_ensemble(outcome_type, check_variation(training[i, ], outcome, learners))

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # predictions
      m_natural$train[jt, tau] <- bound(predict_sl3(fit, nshift_task))
      m_shifted$train[jt, tau] <- bound(predict_sl3(fit, shift_task))
      m_natural$valid[jv, tau] <- bound(predict_sl3(fit, vnshift_task))
      m_shifted$valid[jv, tau] <- bound(predict_sl3(fit, vshift_task))

    }

    # recursion
    estimate_sdr(training = training,
                 shifted = shifted,
                 validation = validation,
                 validation_shifted = validation_shifted,
                 outcome = pseudo,
                 node_list = node_list,
                 C = C,
                 tau = tau - 1,
                 max = max,
                 outcome_type = "continuous",
                 learners = learners,
                 m_shifted = m_shifted,
                 m_natural = m_natural,
                 r = r,
                 pb = pb)
  } else {
    # returns
    out <- list(natural = m_natural$valid,
                shifted = m_shifted$valid)

    return(out)
  }

}
