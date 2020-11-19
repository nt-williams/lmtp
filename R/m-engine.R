estimate_sub <- function(training, shifted, validation, validation_shifted, outcome,
                         node_list, C, deterministic, tau, outcome_type,
                         learners = NULL, m, pb, sl_weights) {
  if (tau > 0) {
    i <- create_censoring_indicators(training, C, tau)$i
    jt <- create_censoring_indicators(training, C, tau)$j
    jv <- create_censoring_indicators(validation, C, tau)$j
    dt <- create_determ_indicators(training, deterministic, tau)
    dv <- create_determ_indicators(validation, deterministic, tau)
    pseudo <- paste0("psi", tau)

    fit <- run_ensemble(training[i & !dt, ][[outcome]],
                        training[i & !dt, node_list[[tau]]],
                        check_variation(training[i & !dt, ][[outcome]],
                                        learners),
                        outcome_type,
                        id = training[i & !dt, ][["lmtp_id"]])
    sl_weights[[tau]] <- extract_sl_weights(fit)

    training[jt & !dt, pseudo] <- bound(predict(fit, shifted[jt & !dt, node_list[[tau]]])$pred)
    training[dt, pseudo] <- 1

    m[jv & !dv, tau] <- bound(predict(fit, validation_shifted[jv & !dv, node_list[[tau]]])$pred)
    m[dv, tau] <- 1

    pb()
    estimate_sub(
      training = training,
      shifted = shifted,
      validation = validation,
      validation_shifted = validation_shifted,
      outcome = pseudo,
      node_list = node_list,
      C = C,
      deterministic = deterministic,
      tau = tau - 1,
      outcome_type = "continuous",
      learners = learners,
      m = m,
      pb = pb,
      sl_weights = sl_weights
    )
  } else {
    list(m = m, sl_weights = sl_weights)
  }
}

estimate_tmle <- function(training, shifted, validation, validation_shifted,
                          outcome, node_list, C, deterministic, tau, outcome_type,
                          m_natural, m_shifted, r, learners = NULL, pb, sl_weights) {
  if (tau > 0) {
    i <- create_censoring_indicators(training, C, tau)$i
    jt <- create_censoring_indicators(training, C, tau)$j
    jv <- create_censoring_indicators(validation, C, tau)$j
    dt <- create_determ_indicators(training, deterministic, tau)
    dv <- create_determ_indicators(validation, deterministic, tau)
    pseudo <- paste0("psi", tau)

    fit <- run_ensemble(training[i & !dt, ][[outcome]],
                        training[i & !dt, node_list[[tau]]],
                        check_variation(training[i & !dt, ][[outcome]],
                                        learners),
                        outcome_type,
                        id = training[i & !dt, ][["lmtp_id"]])
    sl_weights[[tau]] <- extract_sl_weights(fit)

    m_natural$train[jt & !dt, tau] <- bound(predict(fit, training[jt & !dt, node_list[[tau]]])$pred)
    m_shifted$train[jt & !dt, tau] <- bound(predict(fit, shifted[jt & !dt, node_list[[tau]]])$pred)
    m_natural$valid[jv & !dv, tau] <- bound(predict(fit, validation[jv & !dv, node_list[[tau]]])$pred)
    m_shifted$valid[jv & !dv, tau] <- bound(predict(fit, validation_shifted[jv & !dv, node_list[[tau]]])$pred)

    fit <- sw(glm(training[i & !dt, ][[outcome]] ~ offset(qlogis(m_natural$train[i & !dt, tau])),
                                weights = r$train[i & !dt, tau], family = "binomial"))

    training[, pseudo] <- bound(plogis(qlogis(m_shifted$train[, tau]) + coef(fit)))
    training[dt, pseudo] <- 1

    m_natural$valid[, tau] <- bound(plogis(qlogis(m_natural$valid[, tau]) + coef(fit)))
    m_shifted$valid[, tau] <- bound(plogis(qlogis(m_shifted$valid[, tau]) + coef(fit)))
    m_natural$valid[dv, tau] <- 1
    m_shifted$valid[dv, tau] <- 1

    pb()
    estimate_tmle(
      training = training,
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
      sl_weights = sl_weights
    )
  } else {
    list(natural = m_natural$valid,
         shifted = m_shifted$valid,
         sl_weights = sl_weights)
  }
}

estimate_sdr <- function(training, shifted, validation, validation_shifted,
                         outcome, node_list, C, deterministic, tau, max, outcome_type,
                         learners = NULL, m_shifted, m_natural, r, pb, sl_weights) {
  if (tau > 0) {
    i <- create_censoring_indicators(training, C, tau)$i
    jt <- create_censoring_indicators(training, C, tau)$j
    jv <- create_censoring_indicators(validation, C, tau)$j
    dt <- create_determ_indicators(training, deterministic, tau)
    dv <- create_determ_indicators(validation, deterministic, tau)
    pseudo <- paste0("psi", tau + 1)

    if (tau == max) {
      fit <- run_ensemble(training[i & !dt, ][[outcome]],
                          training[i & !dt, node_list[[tau]]],
                          check_variation(training[i & !dt, ][[outcome]],
                                          learners),
                          outcome_type,
                          id = training[i & !dt, ][["lmtp_id"]])
      sl_weights[[tau]] <- extract_sl_weights(fit)

      m_natural$train[jt & !dt, tau] <- bound(predict(fit, training[jt & !dt, node_list[[tau]]])$pred)
      m_shifted$train[jt & !dt, tau] <- bound(predict(fit, shifted[jt & !dt, node_list[[tau]]])$pred)
      m_natural$train[dt, tau] <- 1
      m_shifted$train[dt, tau] <- 1

      m_natural$valid[jv & !dv, tau] <- bound(predict(fit, validation[jv & !dv, node_list[[tau]]])$pred)
      m_shifted$valid[jv & !dv, tau] <- bound(predict(fit, validation_shifted[jv & !dv, node_list[[tau]]])$pred)
      m_natural$valid[dv, tau] <- 1
      m_shifted$valid[dv, tau] <- 1
    }

    if (tau < max) {
      training[, pseudo]  <-
        shifted[, pseudo] <-
        transform_sdr(ratio_sdr(r$train, tau, max),
                      tau, max, m_shifted$train, m_natural$train)

      fit <- run_ensemble(training[i & !dt, ][[pseudo]],
                          training[i & !dt, node_list[[tau]]],
                          check_variation(training[i & !dt, ][[pseudo]],
                                          learners),
                          outcome_type,
                          id = training[i & !dt, ][["lmtp_id"]])
      sl_weights[[tau]] <- extract_sl_weights(fit)

      m_natural$train[jt & !dt, tau] <- bound(predict(fit, training[jt & !dt, node_list[[tau]]])$pred)
      m_shifted$train[jt & !dt, tau] <- bound(predict(fit, shifted[jt & !dt, node_list[[tau]]])$pred)
      m_natural$train[dt, tau] <- 1
      m_shifted$train[dt, tau] <- 1

      m_natural$valid[jv & !dv, tau] <- bound(predict(fit, validation[jv & !dv, node_list[[tau]]])$pred)
      m_shifted$valid[jv & !dv, tau] <- bound(predict(fit, validation_shifted[jv & !dv, node_list[[tau]]])$pred)
      m_natural$valid[dv, tau] <- 1
      m_shifted$valid[dv, tau] <- 1
    }

    pb()
    estimate_sdr(
      training = training,
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
      sl_weights = sl_weights
    )
  } else {
    list(natural = m_natural$valid,
         shifted = m_shifted$valid,
         sl_weights = sl_weights)
  }
}
