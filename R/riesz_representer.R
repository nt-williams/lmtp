cf_rr <- function(task, learners, mtp, control, progress_bar) {
  out <- list()
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- #future::future({
      estimate_rr(get_folded_data(task$natural, task$folds, fold),
                 get_folded_data(task$shifted, task$folds, fold),
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
    #},
    #seed = TRUE)
  }

  out <- future::value(out)
  recombine_ratios(out, task$folds)
}


estimate_rr <- function(natural, shifted, trt, cens, risk, tau, conditional, node_list, learners, mtp, control, progress_bar) {
  representers <- matrix(nrow = nrow(natural$valid), ncol = tau)
  fits <- list()

  for (t in 1:tau) {
    jrt <- censored(natural$train, cens, t)$j
    drt <- at_risk(natural$train, risk, t)
    irv <- censored(natural$valid, cens, t)$i
    jrv <- censored(natural$valid, cens, t)$j
    drv <- at_risk(natural$valid, risk, t)

    trt_t <- ifelse(length(trt) > 1, trt[t], trt)

    frv <- followed_rule(natural$valid[[trt_t]], shifted$valid[[trt_t]], mtp)

    vars <- c(node_list[[t]], cens[[t]])

    fit <- run_riesz_ensemble(
      learners,
      natural$train[jrt & drt, vars, drop = FALSE],
      shifted$train[jrt & drt, vars, drop = FALSE],
      conditional$train[jrt & drt, t, drop = FALSE],
      natural$valid[jrv & drv, vars, drop = FALSE],
      shifted$valid[jrv & drv, vars, drop = FALSE],
      conditional$valid[jrv & drv, t, drop = FALSE],
      folds = control$.learners_trt_folds
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred[jrv & drv, ] <- fit$predictions

    representers[, t] <- pred

    progress_bar()
  }

  list(ratios = representers, fits = fits)
}
