cf_rr <- function(task, G, learners, mtp, control, progress_bar) {
  out <- list()
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_rr(get_folded_data(task$natural, task$folds, fold),
                 get_folded_data(task$shifted, task$folds, fold),
                 get_folded_data(G, task$folds, fold),
                 task$trt,
                 task$cens,
                 task$risk,
                 task$tau,
                 list(train = task$conditional[task$folds[[fold]]$training_set, , drop = FALSE],
                      valid = task$conditional[task$folds[[fold]]$validation_set, , drop = FALSE]),
                 task$node_list$trt,
                 learners,
                 mtp,
                 control,
                 progress_bar)
    },
    seed = TRUE)
  }

  out <- future::value(out)
  recombine_ratios(out, task$folds)
}


estimate_rr <- function(natural, shifted, G, trt, cens, risk, tau, conditional, node_list, learners, mtp, control, progress_bar) {
  representers <- matrix(nrow = nrow(natural$valid), ncol = tau)
  fits <- list()

  prev_riesz <- matrix(1, ncol = 1, nrow = nrow(natural$train))

  for (t in 1:tau) {
    jrt <- censored(natural$train, cens, t)$j
    drt <- at_risk(natural$train, risk, t)
    irv <- censored(natural$valid, cens, t)$i
    jrv <- censored(natural$valid, cens, t)$j
    drv <- at_risk(natural$valid, risk, t)

    if (length(trt) > 1) {
      trt_t <- trt[[t]]
    } else {
      trt_t <- trt[[1]]
    }

    frv <- followed_rule(natural$valid[, trt_t], shifted$valid[, trt_t], mtp)

    vars <- c(node_list[[t]], cens[[t]])

    cumulative_indicator_train <- matrix(as.logical(apply(conditional$train[, (t + 1):(tau + 1), drop = FALSE], 1, prod)), ncol = 1)
    cumulative_indicator_valid <- matrix(as.logical(apply(conditional$valid[, (t + 1):(tau + 1), drop = FALSE], 1, prod)), ncol = 1)

    new_shifted <- natural$train
    new_shifted[[trt_t]] <- shifted$train[[trt_t]]

    new_shifted_valid <- natural$valid
    new_shifted_valid[[trt_t]] <- shifted$valid[[trt_t]]

    fit <- run_riesz_ensemble(
      learners,
      natural$train[jrt & drt, vars, drop = FALSE],
      new_shifted[jrt & drt, vars, drop = FALSE],
      cumulative_indicator_train[jrt & drt, drop = FALSE],
      G$train[jrt & drt, t + 1, drop = FALSE],
      natural$valid[jrv & drv, vars, drop = FALSE],
      new_shifted_valid[jrv & drv, vars, drop = FALSE],
      cumulative_indicator_valid[jrv & drv, drop = FALSE],
      G$valid[jrv & drv, t + 1, drop = FALSE],
      prev_riesz,
      folds = control$.learners_trt_folds
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred[jrv & drv, ] <- fit$predictions
    prev_riesz <- matrix(fit$predictions_train, ncol = 1)

    representers[, t] <- pred

    progress_bar()
  }

  list(ratios = representers, fits = fits)
}
