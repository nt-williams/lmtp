estimate_sub <- function(training, shifted, validation, validation_shifted, outcome,
                         node_list, cens, risk, tau, outcome_type,
                         learners, m, pb, sl_weights, SL_folds) {
  if (tau > 0) {
    i <- censored(training, cens, tau)$i
    jt <- censored(training, cens, tau)$j
    jv <- censored(validation, cens, tau)$j
    rt <- at_risk(training, risk, tau)
    rv <- at_risk(validation, risk, tau)

    pseudo <- paste0("psi", tau)
    vars <- node_list[[tau]]

    fit_task <- initiate_sl3_task(
      training[i & rt, ],
      outcome,
      vars,
      outcome_type,
      "lmtp_id",
      SL_folds
    )

    shift_task <- initiate_sl3_task(
      shifted[jt & rt, ],
      NULL,
      vars,
      NULL,
      "lmtp_id",
      SL_folds
    )

    valid_task <- initiate_sl3_task(
      validation_shifted[jv & rv, ],
      NULL,
      vars,
      NULL,
      "lmtp_id",
      SL_folds
    )

    ensemble <- initiate_ensemble(
      outcome_type,
      check_variation(training[i & rt, ][[outcome]], learners)
    )

    fit <- run_ensemble(ensemble, fit_task)

    sl_weights[[tau]] <- extract_sl_weights(fit)

    training[jt & rt, pseudo] <- SL_predict(fit, shift_task)
    m[jv & rv, tau] <- SL_predict(fit, valid_task)

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
      cens = cens,
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
                          outcome, node_list, cens, risk, tau, outcome_type,
                          m_natural, m_shifted, ratios, learners, pb,
                          weights, sl_weights, SL_folds) {
  if (tau > 0) {
    i <- censored(training, cens, tau)$i
    jt <- censored(training, cens, tau)$j
    jv <- censored(validation, cens, tau)$j

    rt <- at_risk(training, risk, tau)
    rv <- at_risk(validation, risk, tau)

    pseudo <- paste0("psi", tau)
    vars <- node_list[[tau]]

    fit_task <-
      initiate_sl3_task(training[i & rt,], outcome, vars, outcome_type, "lmtp_id", SL_folds)

    nshift_task <-
      initiate_sl3_task(training[jt & rt,], NULL, vars, NULL, "lmtp_id", SL_folds)

    shift_task <-
      initiate_sl3_task(shifted[jt & rt,], NULL, vars, NULL, "lmtp_id", SL_folds)

    vnshift_task <-
      initiate_sl3_task(validation[jv & rv,], NULL, vars, NULL, "lmtp_id", SL_folds)

    vshift_task <-
      initiate_sl3_task(validation_shifted[jv & rv,], NULL, vars, NULL, "lmtp_id", SL_folds)

    ensemble <-
      initiate_ensemble(
        outcome_type,
        check_variation(training[i & rt,][[outcome]], learners)
      )

    fit <- run_ensemble(ensemble, fit_task)

    sl_weights[[tau]] <- extract_sl_weights(fit)

    m_natural$train[jt & rt, tau] <- SL_predict(fit, nshift_task)
    m_shifted$train[jt & rt, tau] <- SL_predict(fit, shift_task)
    m_natural$valid[jv & rv, tau] <- SL_predict(fit, vnshift_task)
    m_shifted$valid[jv & rv, tau] <- SL_predict(fit, vshift_task)

    wts <- {
      if (is.null(weights))
        ratios[i & rt, tau]
      else
        ratios[i & rt, tau] * weights[i & rt]
    }

    fit <- sw(
      glm(
        training[i & rt, ][[outcome]]
        ~ offset(qlogis(m_natural$train[i & rt, tau])),
        weights = wts,
        family = "binomial"
      )
    )

    training[jt & rt, pseudo] <-
      bound(
        plogis(qlogis(m_shifted$train[jt & rt, tau]) + coef(fit))
      )

    m_natural$valid[jv & rv, tau] <-
      bound(
        plogis(qlogis(m_natural$valid[jv & rv, tau]) + coef(fit))
      )

    m_shifted$valid[jv & rv, tau] <-
      bound(
        plogis(qlogis(m_shifted$valid[jv & rv, tau]) + coef(fit))
      )

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
      cens = cens,
      risk = risk,
      tau = tau - 1,
      outcome_type = "continuous",
      m_natural = m_natural,
      m_shifted = m_shifted,
      ratios = ratios,
      learners = learners,
      pb = pb,
      weights = weights,
      sl_weights = sl_weights,
      SL_folds = SL_folds
    )
  } else {
    list(
      natural = m_natural$valid,
      shifted = m_shifted$valid,
      sl_weights = sl_weights
    )
  }
}

estimate_sdr <- function(training, shifted, validation, validation_shifted,
                         outcome, node_list, cens, risk, tau, max, outcome_type,
                         learners, m_shifted, m_natural, ratios,
                         pb, sl_weights, SL_folds) {
  if (tau > 0) {
    i <- censored(training, cens, tau)$i
    jt <- censored(training, cens, tau)$j
    jv <- censored(validation, cens, tau)$j
    rt <- at_risk(training, risk, tau)
    rv <- at_risk(validation, risk, tau)

    pseudo <- paste0("psi", tau + 1)
    vars <- node_list[[tau]]

    if (tau == max) {
      fit_task <- initiate_sl3_task(
        training[i & rt, ],
        outcome,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      nshift_task <- initiate_sl3_task(
        training[jt & rt, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      shift_task <- initiate_sl3_task(
        shifted[jt & rt, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      vnshift_task <- initiate_sl3_task(
        validation[jv & rv, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      vshift_task <- initiate_sl3_task(
        validation_shifted[jv & rv, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      ensemble <- initiate_ensemble(
        outcome_type,
        check_variation(training[i & rt, ][[outcome]], learners)
      )

      fit <- run_ensemble(ensemble, fit_task)
      sl_weights[[tau]] <- extract_sl_weights(fit)

      m_natural$train[jt & rt, tau] <- SL_predict(fit, nshift_task)
      m_shifted$train[jt & rt, tau] <- SL_predict(fit, shift_task)
      m_natural$valid[jv & rv, tau] <- SL_predict(fit, vnshift_task)
      m_shifted$valid[jv & rv, tau] <- SL_predict(fit, vshift_task)

      m_natural$train[!rt, tau] <- 0
      m_shifted$train[!rt, tau] <- 0
      m_natural$valid[!rv, tau] <- 0
      m_shifted$valid[!rv, tau] <- 0
    }

    if (tau < max) {
      training[, pseudo]  <-
        shifted[, pseudo] <-
        transform_sdr(
          ratio_sdr(ratios, tau, max),
          tau, max,
          m_shifted$train,
          m_natural$train
        )

      fit_task <- initiate_sl3_task(
        training[i & rt, ],
        pseudo,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      nshift_task <- initiate_sl3_task(
        training[jt & rt, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      shift_task <- initiate_sl3_task(
        shifted[jt & rt, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      vnshift_task <- initiate_sl3_task(
        validation[jv & rv, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      vshift_task <- initiate_sl3_task(
        validation_shifted[jv & rv, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      ensemble <- initiate_ensemble(
        outcome_type,
        check_variation(training[i & rt, ][[pseudo]], learners)
      )

      fit <- run_ensemble(ensemble, fit_task)

      sl_weights[[tau]] <- extract_sl_weights(fit)

      m_natural$train[jt & rt, tau] <- SL_predict(fit, nshift_task)
      m_shifted$train[jt & rt, tau] <- SL_predict(fit, shift_task)
      m_natural$valid[jv & rv, tau] <- SL_predict(fit, vnshift_task)
      m_shifted$valid[jv & rv, tau] <- SL_predict(fit, vshift_task)

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
      cens = cens,
      risk = risk,
      tau = tau - 1,
      max = max,
      outcome_type = "continuous",
      learners = learners,
      m_shifted = m_shifted,
      m_natural = m_natural,
      ratios = ratios,
      pb = pb,
      sl_weights = sl_weights,
      SL_folds = SL_folds
    )
  } else {
    list(
      natural = m_natural$valid,
      shifted = m_shifted$valid,
      sl_weights = sl_weights
    )
  }
}
