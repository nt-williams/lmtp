estimate_r <- function(data, shifted, trt, cens, risk, tau,
                       node_list, learners, pb, sl_weights,
                       intervention_type, SL_folds) {
  training <- data$train
  validation <- data$valid

  nt <- nrow(training)
  nv <- nrow(validation)

  returns <- matrix(nrow = nv, ncol = tau)

  for (t in 1:tau) {
    jrt <- rep(censored(training, cens, t)$j, 2)
    drt <- rep(at_risk(training, risk, t), 2)
    irv <- rep(censored(validation, cens, t)$i, 2)
    jrv <- rep(censored(validation, cens, t)$j, 2)
    drv <- rep(at_risk(validation, risk, t), 2)

    frv <- followed_rule(
      validation[[trt[t]]],
      shifted$valid[[trt[t]]],
      intervention_type
    )

    vars <- c(node_list[[t]], cens[[t]])
    stcks <- create_r_stacks(data, shifted, trt, cens, t)

    fit <- run_ensemble(
      stcks$train[jrt & drt, ]$si,
      stcks$train[jrt & drt, vars],
      learners,
      "binomial",
      stcks$train[jrt & drt, ]$lmtp_id,
      SL_folds
    )

    sl_weights[[t]] <- extract_sl_weights(fit)

    pred <- matrix(-999L, nrow = nv * 2, ncol = 1)

    pred[jrv & drv, ] <- SL_predict(
      fit, stcks$valid[jrv & drv, vars],
      .Machine$double.eps
    )

    ratios <- create_ratios(
      pred, irv, drv, frv,
      intervention_type == "mtp"
    )

    returns[, t] <- ratios[stcks$valid$si == 0]

    pb()
  }

  list(ratios = returns, sl_weights = sl_weights)
}

create_ratios <- function(pred, cens, risk, followed, mtp) {
  use <- ifelse(followed & isFALSE(mtp), pmax(pred, 0.5), pred)
  (use * cens * risk * followed) / (1 - pmin(use, 0.99))
}

ratio_tmle_ipw <- function(ratios) {
  matrix(
    t(apply(ratios[["ratios"]], 1, cumprod)),
    nrow = nrow(ratios[["ratios"]]),
    ncol = ncol(ratios[["ratios"]])
  )
}

ratio_sdr <- function(ratio, tau, max_tau) {
  out <- t(
    apply(
      ratio[, (tau + 1):max_tau, drop = FALSE],
      1,
      cumprod
    )
  )

  if (tau != max_tau - 1) {
    return(out)
  }

  t(out)
}
