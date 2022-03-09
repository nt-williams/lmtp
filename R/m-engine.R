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
      fit_task <- sl3_task(
        training[i & rt, ],
        outcome,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      nshift_task <- sl3_task(
        training[jt & rt, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      shift_task <- sl3_task(
        shifted[jt & rt, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      vnshift_task <- sl3_task(
        validation[jv & rv, ],
        NULL,
        vars,
        NULL,
        "lmtp_id",
        SL_folds
      )

      vshift_task <- sl3_task(
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

      fit_task <- sl3_task(
        training[i & rt, ],
        pseudo,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      nshift_task <- sl3_task(
        training[jt & rt, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      shift_task <- sl3_task(
        shifted[jt & rt, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      vnshift_task <- sl3_task(
        validation[jv & rv, ],
        NULL,
        vars,
        outcome_type,
        "lmtp_id",
        SL_folds
      )

      vshift_task <- sl3_task(
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
