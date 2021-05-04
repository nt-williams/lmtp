estimate_r <- function(data, shifted, trt, cens, risk, tau,
                       node_list, learners, pb, sl_weights,
                       intervention_type, SL_folds) {
  training <- data$train
  validation <- data$valid
  nt <- nrow(training)
  nv <- nrow(validation)
  rt <- list(natural = matrix(nrow = nt, ncol = tau),
             shifted = matrix(nrow = nt, ncol = tau))
  rv <- list(natural = matrix(nrow = nv, ncol = tau),
             shifted = matrix(nrow = nv, ncol = tau))

  for (t in 1:tau) {
    irt <- rep(censored(training, cens, t)$i, 2)
    jrt <- rep(censored(training, cens, t)$j, 2)
    drt <- rep(at_risk(training, risk, t), 2)
    irv <- rep(censored(validation, cens, t)$i, 2)
    jrv <- rep(censored(validation, cens, t)$j, 2)
    drv <- rep(at_risk(validation, risk, t), 2)
    frt <- followed_rule(training[[trt[t]]], shifted$train[[trt[t]]], intervention_type)
    frv <- followed_rule(validation[[trt[t]]], shifted$valid[[trt[t]]], intervention_type)
    vars <- c(node_list[[t]], cens[[t]])
    stcks <- create_r_stacks(data, shifted, trt, cens, t)

    fit_task <- initiate_sl3_task(
      stcks$train[jrt & drt,], "si",
      vars, "binomial", "lmtp_id", SL_folds
    )
    vpred_task <- initiate_sl3_task(
      stcks$valid[jrv & drv, ], "si",
      vars, "binomial", "lmtp_id", SL_folds
    )

    ensemble <- initiate_ensemble("binomial", learners)
    fit <- run_ensemble(ensemble, fit_task)
    sl_weights[[t]] <- extract_sl_weights(fit)

    pred <- matrix(-999L, nrow = nt * 2, ncol = 1)
    pred[jrt & drt, ] <- SL_predict(fit, fit_task, .Machine$double.eps)

    rat <- create_ratios(pred, irt, drt, frt, intervention_type == "mtp")
    rt$natural[, t] <- rat[stcks$train$si == 0]
    rt$shifted[, t] <- rat[stcks$train$si == 1]

    pred <- matrix(-999L, nrow = nv * 2, ncol = 1)
    pred[jrv & drv, ] <- SL_predict(fit, vpred_task, .Machine$double.eps)

    rat <- create_ratios(pred, irv, drv, frv, intervention_type == "mtp")
    rv$natural[, t] <- rat[stcks$valid$si == 0]
    rv$shifted[, t] <- rat[stcks$valid$si == 1]

    pb()
  }
  list(train = rt,
       valid = rv,
       sl_weights = sl_weights)
}

create_ratios <- function(pred, cens, risk, followed, mtp) {
  use <- ifelse(followed & isFALSE(mtp), pmax(pred, 0.5), pred)
  (use * cens * risk * followed) / (1 - use)
}

ratio_dr <- function(ratios, V, trim) {
  out <- list()
  for (i in 1:V) {
    out[[i]] <- list()
    out[[i]]$train <- check_extreme_ratio(
      matrix(t(apply(ratios[[i]]$train$natural, 1, cumprod)),
             nrow = nrow(ratios[[i]]$train$natural),
             ncol = ncol(ratios[[i]]$train$natural)),
      trim
    )
    out[[i]]$valid <- check_extreme_ratio(
      matrix(t(apply(ratios[[i]]$valid$natural, 1, cumprod)),
             nrow = nrow(ratios[[i]]$valid$natural),
             ncol = ncol(ratios[[i]]$valid$natural)),
      trim
    )
    out[[i]]$sl_weights <- ratios[[i]]$sl_weights
  }
  return(out)
}

ratio_ipw <- function(ratio, trim) {
    list(r = check_extreme_ratio(
      matrix(
        t(apply(ratio$r, 1, cumprod)),
        nrow = nrow(ratio$r),
        ncol = ncol(ratio$r)
      ),
      trim
    ),
    sl_weights = ratio$sl_weights)
}

ratio_sdr <- function(ratio, tau, max_tau, trim) {
  out <- t(apply(ratio$natural[, (tau + 1):max_tau, drop = FALSE], 1, cumprod))
  if (tau == max_tau - 1) out <- t(out)
  check_extreme_ratio(out, trim)
}
