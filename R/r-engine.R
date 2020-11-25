estimate_r <- function(training, validation, trt, cens, risk, shift,
                       tau, node_list, learners = NULL, pb, sl_weights) {
  nt <- nrow(training)
  nv <- nrow(validation)
  rt <- list(natural = matrix(nrow = nt, ncol = tau),
             shifted = matrix(nrow = nt, ncol = tau))
  rv <- list(natural = matrix(nrow = nv, ncol = tau),
             shifted = matrix(nrow = nv, ncol = tau))

  for (t in 1:tau) {
    irt <- rep(censored(training, cens, t)$j, 2)
    drt <- rep(at_risk(training, risk, t), 2)
    irv <- rep(censored(validation, cens, t)$j, 2)
    drv <- rep(at_risk(validation, risk, t), 2)

    stcks <- create_r_stacks(training, validation, trt, cens, shift, t, nt, nv)
    fit <- run_ensemble(subset(stcks$train, irt & drt)$si,
                        subset(stcks$train, irt & drt)[, c(node_list[[t]], cens[[t]])],
                        learners,
                        "binomial",
                        subset(stcks$train, irt & drt)$lmtp_id)
    sl_weights[[t]] <- extract_sl_weights(fit)

    pred <- matrix(nrow = nt * 2, ncol = 1)
    pred[irt & drt, ] <- bound(
      predict(fit, stcks$train[irt & drt, c(node_list[[t]], cens[[t]])])$pred,
      .Machine$double.eps
    )

    rat <- create_ratios(pred, training, cens, t)
    rt$natural[, t] <- rat[stcks$train$si == 0]
    rt$shifted[, t] <- rat[stcks$train$si == 1]

    pred <- matrix(nrow = nv * 2, ncol = 1)
    pred[irv & drv, ] <- bound(
      predict(fit, stcks$valid[irv & drv, c(node_list[[t]], cens[[t]])])$pred,
      .Machine$double.eps
    )

    rat <- create_ratios(pred, validation, cens, t)
    rv$natural[, t] <- rat[stcks$valid$si == 0]
    rv$shifted[, t] <- rat[stcks$valid$si == 1]

    pb()
  }
  list(train = rt,
       valid = rv,
       sl_weights = sl_weights)
}

create_ratios <- function(pred, data, cens, tau) {
  out <- pred * rep(censored(data, cens, tau)$i, 2) / (1 - pred)
  out <- ifelse(is.na(out), -999, out)
  return(out)
}

# create_ratios <- function(pred) {
#   out <- pred / (1 - pred)
#   ifelse(is.na(out), 0, out)
# }

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
