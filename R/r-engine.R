
estimate_r <- function(training, validation, trt, cens, deterministic, shift,
                       tau, node_list, learners = NULL, pb, sl_weights) {
  nt <- nrow(training)
  nv <- nrow(validation)
  rt <- list(natural = matrix(nrow = nt, ncol = tau),
             shifted = matrix(nrow = nt, ncol = tau))
  rv <- list(natural = matrix(nrow = nv, ncol = tau),
             shifted = matrix(nrow = nv, ncol = tau))

  for (t in 1:tau) {
    i     <- rep(create_censoring_indicators(training, cens, t)$j, 2)
    d     <- rep(create_determ_indicators(training, deterministic, t), 2)
    stcks <- create_r_stacks(training, validation, trt, cens, shift, t, nt, nv)

    fit <- run_ensemble(subset(stcks$train, i & !d)$si,
                        subset(stcks$train, i & !d)[, c(node_list[[t]], cens[[t]])],
                        learners,
                        "binomial",
                        subset(stcks$train, i & !d)$lmtp_id)
    sl_weights[[t]] <- extract_sl_weights(fit)

    pred <- bound(predict(fit, stcks$train[, c(node_list[[t]], cens[[t]])])$pred,
                  .Machine$double.eps)
    rat <- create_ratios(pred, training, cens, t)
    rt$natural[, t] <- rat[stcks$train$si == 0]
    rt$shifted[, t] <- rat[stcks$train$si == 1]
    rt$natural[create_determ_indicators(training, deterministic, t), t] <- 1
    rt$shifted[create_determ_indicators(training, deterministic, t), t] <- 1

    pred <- bound(predict(fit, stcks$valid[, c(node_list[[t]], cens[[t]])])$pred,
                  .Machine$double.eps)
    rat <- create_ratios(pred, validation, cens, t)
    rv$natural[, t] <- rat[stcks$valid$si == 0]
    rv$shifted[, t] <- rat[stcks$valid$si == 1]
    rv$natural[create_determ_indicators(validation, deterministic, t), t] <- 1
    rv$shifted[create_determ_indicators(validation, deterministic, t), t] <- 1

    pb()
  }
  list(train = rt,
       valid = rv,
       sl_weights = sl_weights)
}

create_ratios <- function(pred, data, cens, tau) {
  out <- pred * rep(create_censoring_indicators(data, cens, tau)$i, 2) / (1 - bound(pred))
  out <- ifelse(is.na(out), 0, out)
  return(out)
}

ratio_dr <- function(ratios, V) {
  out <- list()
  for (i in 1:V) {
      out[[i]] <- list()
      out[[i]]$train <- check_extreme_ratio(
        matrix(t(apply(ratios[[i]]$train$natural, 1, cumprod)),
               nrow = nrow(ratios[[i]]$train$natural),
               ncol = ncol(ratios[[i]]$train$natural))
      )
      out[[i]]$valid <- check_extreme_ratio(
        matrix(t(apply(ratios[[i]]$valid$natural, 1, cumprod)),
               nrow = nrow(ratios[[i]]$valid$natural),
               ncol = ncol(ratios[[i]]$valid$natural))
      )
      out[[i]]$sl_weights <- ratios[[i]]$sl_weights
  }
  return(out)
}

ratio_ipw <- function(ratio) {
  out <-
    list(r = check_extreme_ratio(matrix(
      t(apply(ratio$r, 1, cumprod)),
      nrow = nrow(ratio$r),
      ncol = ncol(ratio$r)
    )),
    sl_weights = ratio$sl_weights)
  return(out)
}

ratio_sdr <- function(ratio, tau, max_tau) {
  out <- t(apply(ratio$natural[, (tau + 1):max_tau, drop = FALSE], 1, cumprod))
  if (tau == max_tau - 1) out <- t(out)
  return(check_extreme_ratio(out))
}
