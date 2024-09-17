cf_G <- function(task, learners, control, pb) {
  out <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_G(
        get_folded_data(task$natural, task$folds, fold),
        get_folded_data(task$shifted, task$folds, fold),
        get_folded_data(task$conditional, task$folds, fold),
        task$trt,
        task$cens,
        task$risk,
        task$tau,
        task$node_list$trt,
        learners,
        pb,
        mtp,
        control
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    G = recombine_outcome(out, "G", task$folds),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_G <- function(natural,
                       shifted,
                       conditional,
                       trt,
                       cens,
                       risk,
                       tau,
                       node_list,
                       learners,
                       pb,
                       mtp,
                       control) {
  G <- matrix(nrow = nrow(natural$valid), ncol = tau + 1)
  G[, tau + 1] <- 1
  fits <- vector("list", length = tau)

  for (t in 0:(tau - 1)) {
    if (t == 0) {
      cumulative_indicator <- apply(conditional$train[, 1:(tau + 1), drop = FALSE], 1, prod)
      G[, t + 1] <- mean(cumulative_indicator)
      next
    }

    jrt <- censored(natural$train, cens, t)$j
    jrv <- censored(natural$valid, cens, t)$j
    drt <- at_risk(natural$train, risk, t)
    drv <- at_risk(natural$valid, risk, t)

    vars <- c(node_list[[t]], cens[[t]])

    train <- natural$train[jrt & drt, c("lmtp_id", vars)]
    train$tmp_lmtp_pseudo <- apply(conditional$train[, (t + 1):(tau + 1), drop = FALSE], 1, prod)

    valid <- natural$valid[jrv & drv, c("lmtp_id", vars)]

    if (all(train$tmp_lmtp_pseudo == 1)) {
      G[, t + 1] <- 1
      next
    }

    if (all(data_train$tmp_lmtp_pseudo == 0)) {
      G[, t + 1] <- 0
      next
    }

    fit <- run_ensemble(
      train,
      "tmp_lmtp_pseudo",
      learners,
      "binomial",
      "lmtp_id",
      control$.learners_conditional_folds
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred <- bound(SL_predict(fit, data_valid[, c("lmtp_id", vars)]), .Machine$double.eps)

    G[, t + 1] <- pred

    pb()
  }

  list(G = G, fits = fits)
}
