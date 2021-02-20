estimate_sub <- function(training, shifted, validation, validation_shifted, outcome,
                         node_list, C, risk, tau, outcome_type,
                         learners = NULL, m, pb, sl_weights, SL_folds) {
  if (tau > 0) {
    i <- censored(training, C, tau)$i
    jt <- censored(training, C, tau)$j
    jv <- censored(validation, C, tau)$j
    rt <- at_risk(training, risk, tau)
    rv <- at_risk(validation, risk, tau)
    pseudo <- paste0("psi", tau)
    vars <- node_list[[tau]]

    fit <- run_ensemble(training[i & rt, ][[outcome]],
                        training[i & rt, vars],
                        check_variation(training[i & rt, ][[outcome]],
                                        learners),
                        outcome_type,
                        id = training[i & rt, ][["lmtp_id"]],
                        SL_folds)
    sl_weights[[tau]] <- extract_sl_weights(fit)

    training[jt & rt, pseudo] <- SL_predict(fit, shifted[jt & rt, vars])
    m[jv & rv, tau] <- SL_predict(fit, validation_shifted[jv & rv, vars])

    training[!rt, pseudo] <- 0
    m[!rv, tau] <- 0

    pb()
    estimate_sub(
      training = training,
      shifted = shifted,
      validation = validation,
      validation_shifted = validation_shifted,
      outcome = pseudo,
      node_list = node_list,
      C = C,
      risk = risk,
      tau = tau - 1,
      outcome_type = "continuous",
      learners = learners,
      m = m,
      pb = pb,
      sl_weights = sl_weights,
      SL_folds = SL_folds
    )
  } else {
    list(m = m, sl_weights = sl_weights)
  }
}

estimate_tmle <- function(training, shifted, validation, validation_shifted,
                          outcome, node_list, C, risk, tau, outcome_type,
                          m_natural, m_shifted, r, learners = NULL, pb,
                          weights = NULL, sl_weights, SL_folds) {
  if (tau > 0) {
    i <- censored(training, C, tau)$i
    jt <- censored(training, C, tau)$j
    jv <- censored(validation, C, tau)$j
    rt <- at_risk(training, risk, tau)
    rv <- at_risk(validation, risk, tau)
    pseudo <- paste0("psi", tau)
    vars <- node_list[[tau]]

    fit <- run_ensemble(training[i & rt, ][[outcome]],
                        training[i & rt, vars],
                        check_variation(training[i & rt, ][[outcome]],
                                        learners),
                        outcome_type,
                        id = training[i & rt, ][["lmtp_id"]],
                        SL_folds)
    sl_weights[[tau]] <- extract_sl_weights(fit)

    m_natural$train[jt & rt, tau] <- SL_predict(fit, training[jt & rt, vars])
    m_shifted$train[jt & rt, tau] <- SL_predict(fit, shifted[jt & rt, vars])
    m_natural$valid[jv & rv, tau] <- SL_predict(fit, validation[jv & rv, vars])
    m_shifted$valid[jv & rv, tau] <- SL_predict(fit, validation_shifted[jv & rv, vars])

    wts <- {
      if (is.null(weights))
        r$train[i & rt, tau]
      else
        r$train[i & rt, tau] * weights[i & rt]
    }

    fit <- sw(glm(training[i & rt, ][[outcome]]
                  ~ offset(qlogis(m_natural$train[i & rt, tau])),
                  weights = wts, family = "binomial"))

    training[jt & rt, pseudo] <- bound(plogis(qlogis(m_shifted$train[jt & rt, tau]) + coef(fit)))
    m_natural$valid[jv & rv, tau] <- bound(plogis(qlogis(m_natural$valid[jv & rv, tau]) + coef(fit)))
    m_shifted$valid[jv & rv, tau] <- bound(plogis(qlogis(m_shifted$valid[jv & rv, tau]) + coef(fit)))

    training[!rt, pseudo] <- 0
    m_natural$valid[!rv, tau] <- 0
    m_shifted$valid[!rv, tau] <- 0

    pb()
    estimate_tmle(
      training = training,
      shifted = shifted,
      validation = validation,
      validation_shifted = validation_shifted,
      outcome = pseudo,
      node_list = node_list,
      C = C,
      risk = risk,
      tau = tau - 1,
      outcome_type = "continuous",
      m_natural = m_natural,
      m_shifted = m_shifted,
      r = r,
      learners = learners,
      pb = pb,
      sl_weights = sl_weights,
      SL_folds = SL_folds
    )
  } else {
    list(natural = m_natural$valid,
         shifted = m_shifted$valid,
         sl_weights = sl_weights)
  }
}

estimate_sdr <- function(training, shifted, validation, validation_shifted,
                         outcome, node_list, C, risk, tau, max, outcome_type,
                         learners = NULL, m_shifted, m_natural, r, pb, sl_weights,
                         trim, SL_folds) {
  if (tau > 0) {
    i <- censored(training, C, tau)$i
    jt <- censored(training, C, tau)$j
    jv <- censored(validation, C, tau)$j
    rt <- at_risk(training, risk, tau)
    rv <- at_risk(validation, risk, tau)
    pseudo <- paste0("psi", tau + 1)
    vars <- node_list[[tau]]

    if (tau == max) {
      fit <- run_ensemble(training[i & rt, ][[outcome]],
                          training[i & rt, vars],
                          check_variation(training[i & rt, ][[outcome]],
                                          learners),
                          outcome_type,
                          id = training[i & rt, ][["lmtp_id"]],
                          SL_folds)
      sl_weights[[tau]] <- extract_sl_weights(fit)

      m_natural$train[jt & rt, tau] <- SL_predict(fit, training[jt & rt, vars])
      m_shifted$train[jt & rt, tau] <- SL_predict(fit, shifted[jt & rt, vars])
      m_natural$valid[jv & rv, tau] <- SL_predict(fit, validation[jv & rv, vars])
      m_shifted$valid[jv & rv, tau] <- SL_predict(fit, validation_shifted[jv & rv, vars])

      m_natural$train[!rt, tau] <- 0
      m_shifted$train[!rt, tau] <- 0
      m_natural$valid[!rv, tau] <- 0
      m_shifted$valid[!rv, tau] <- 0
    }

    if (tau < max) {
      training[, pseudo]  <-
        shifted[, pseudo] <-
        transform_sdr(ratio_sdr(r$train, tau, max, trim),
                      tau, max, m_shifted$train, m_natural$train)

      fit <- run_ensemble(training[i & rt, ][[pseudo]],
                          training[i & rt, vars],
                          check_variation(training[i & rt, ][[pseudo]],
                                          learners),
                          outcome_type,
                          id = training[i & rt, ][["lmtp_id"]],
                          SL_folds)
      sl_weights[[tau]] <- extract_sl_weights(fit)

      m_natural$train[jt & rt, tau] <- SL_predict(fit, training[jt & rt, vars])
      m_shifted$train[jt & rt, tau] <- SL_predict(fit, shifted[jt & rt, vars])
      m_natural$valid[jv & rv, tau] <- SL_predict(fit, validation[jv & rv, vars])
      m_shifted$valid[jv & rv, tau] <- SL_predict(fit, validation_shifted[jv & rv, vars])

      m_natural$train[!rt, tau] <- 0
      m_shifted$train[!rt, tau] <- 0
      m_natural$valid[!rv, tau] <- 0
      m_shifted$valid[!rv, tau] <- 0
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
      risk = risk,
      tau = tau - 1,
      max = max,
      outcome_type = "continuous",
      learners = learners,
      m_shifted = m_shifted,
      m_natural = m_natural,
      r = r,
      pb = pb,
      sl_weights = sl_weights,
      trim = trim,
      SL_folds = SL_folds
    )
  } else {
    list(natural = m_natural$valid,
         shifted = m_shifted$valid,
         sl_weights = sl_weights)
  }
}
