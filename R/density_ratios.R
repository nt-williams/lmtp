cf_r <- function(task, learners, mtp, control, pb) {
  out <- vector("list", length = length(task$folds))
  
  if (length(learners) == 1 && learners == "SL.mean") {
    warning("Using 'SL.mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('SL.glm', 'SL.glmnet')`.",
            call. = FALSE)
  }
  
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_r(
        get_folded_data(task$natural, task$folds, fold),
        get_folded_data(task$shifted, task$folds, fold),
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

  trim_ratios(recombine_ratios(future::value(out), task$folds), control$.trim)
}

estimate_r <- function(natural, shifted, trt, cens, risk, tau, node_list, learners, pb, mtp, control) {
  densratios <- matrix(nrow = nrow(natural$valid), ncol = tau)
  fits <- vector("list", length = tau)

  for (t in 1:tau) {
    jrt <- rep(censored(natural$train, cens, t)$j, 2)
    drt <- rep(at_risk(natural$train, risk, t), 2)
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
    stacked <- stack_data(natural$train, shifted$train, trt, cens, t)

    fit <- run_ensemble(
      stacked[jrt & drt, ][["tmp_lmtp_stack_indicator"]],
      stacked[jrt & drt, vars],
      learners,
      "binomial",
      stacked[jrt & drt, ]$lmtp_id,
      control$.learners_trt_folds
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred[jrv & drv, ] <- bound(SL_predict(fit, natural$valid[jrv & drv, vars]), .Machine$double.eps)

    ratios <- density_ratios(pred, irv, drv, frv, mtp)
    densratios[, t] <- ratios

    pb()
  }

  list(ratios = densratios, fits = fits)
}

stack_data <- function(natural, shifted, trt, cens, tau) {
  shifted_half <- natural

  if (length(trt) > 1 || tau == 1) {
    shifted_half[, trt[[tau]]] <- shifted[, trt[[tau]]]
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
