estimate_r <- function(data, shifted, trt, cens, risk, tau,
                       node_list, learners, pb, sl_weights,
                       intervention_type, SL_folds) {
  training <- data$train
  validation <- data$valid

  nt <- nrow(training)
  nv <- nrow(validation)

  returns <- list(
    natural = matrix(nrow = nv, ncol = tau),
    shifted = matrix(nrow = nv, ncol = tau)
  )

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

    returns$natural[, t] <- ratios[stcks$valid$si == 0]
    returns$shifted[, t] <- ratios[stcks$valid$si == 1]

    pb()
  }

  list(ratios = returns, sl_weights = sl_weights)
}

create_ratios <- function(pred, cens, risk, followed, mtp) {
  use <- ifelse(followed & isFALSE(mtp), pmax(pred, 0.5), pred)
  (use * cens * risk * followed) / (1 - pmin(use, 0.99))
}

ratio_tmle <- function(ratios) {
  matrix(
    t(apply(ratios[["natural"]], 1, cumprod)),
    nrow = nrow(ratios[["natural"]]),
    ncol = ncol(ratios[["natural"]])
  )
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
