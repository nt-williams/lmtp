cf_r <- function(Task, learners, type, lrnr_folds, trim, pb) {
  fopts <- options("lmtp.bound", "lmtp.trt.length")
  out <- list()

  for (fold in seq_along(Task$folds)) {
    out[[fold]] <- future::future({
      options(fopts)

      estimate_r(
        get_folded_data(Task$natural, Task$folds, fold),
        get_folded_data(Task$shifted, Task$folds, fold),
        Task$trt, Task$cens, Task$risk, Task$tau, Task$node_list$trt,
        learners, pb, type, lrnr_folds
      )
    },
    seed = TRUE)
  }

  trim_ratios(recombine_ratios(future::value(out), Task$folds), trim)
}

estimate_r <- function(natural, shifted, trt, cens, risk, tau, node_list, learners, pb, type, lrnr_folds) {
  densratios <- matrix(nrow = nrow(natural$valid), ncol = tau)
  fits <- list()

  for (t in 1:tau) {
    jrt <- rep(censored(natural$train, cens, t)$j, 2)
    drt <- rep(at_risk(natural$train, risk, t), 2)
    irv <- censored(natural$valid, cens, t)$i
    jrv <- censored(natural$valid, cens, t)$j
    drv <- at_risk(natural$valid, risk, t)

    frv <- followed_rule(natural$valid[[trt[t]]], shifted$valid[[trt[t]]], type)

    vars <- c(node_list[[t]], cens[[t]])
    stacked <- stack_data(natural$train, shifted$train, trt, cens, t)

    train_task <- sl3_task(stacked[jrt & drt, ], "tmp_lmtp_stack_indicator", vars, "binomial", "lmtp_id", lrnr_folds)
    pred_task <- sl3_task(natural$valid[jrv & drv, ], NULL, vars, NULL, "lmtp_id")

    fit <- learners$train(train_task)
    fits[[t]] <- fit$fit_object

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred[jrv & drv, ] <- bound(fit$predict(pred_task), .Machine$double.eps)

    ratios <- density_ratios(pred, irv, drv, frv, type == "mtp")
    densratios[, t] <- ratios

    pb()
  }

  list(ratios = densratios, fits = fits)
}

stack_data <- function(natural, shifted, trt, cens, tau) {
  shifted_half <- natural

  if (getOption("lmtp.trt.length") == "standard" || tau == 1) {
    shifted_half[[trt[tau]]] <- shifted[[trt[tau]]]
  }

  if (!is.null(cens)) {
    shifted_half[[cens[tau]]] <- shifted[[cens[tau]]]
  }

  out <- rbind(natural, shifted_half)
  out[["tmp_lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(natural))
  out
}

density_ratios <- function(pred, cens, risk, followed, mtp) {
  pred <- ifelse(followed & isFALSE(mtp), pmax(pred, 0.5), pred)
  (pred * cens * risk * followed) / (1 - pmin(pred, 0.999))
}

ratio_tmle_ipw <- function(ratios) {
  matrix(
    t(apply(ratios, 1, cumprod)),
    nrow = nrow(ratios),
    ncol = ncol(ratios)
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
