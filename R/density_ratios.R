cf_density_ratios <- function(task, learners, mtp, control, pb) {
  ans <- vector("list", length = length(task$folds))

  if (length(learners) == 1 && learners == "mean") {
    warning("Using 'mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('glm', 'cv_glmnet')`.",
            call. = FALSE)
  }

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_density_ratios(task, fold, learners, mtp, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  ans <- list(density_ratios = recombine(rbind_depth(ans, "ratios"), task$folds),
              fits = lapply(ans, function(x) x[["fits"]]))

  ans$density_ratios <- trim(ans$density_ratios, control$.trim)
  ans
}

estimate_density_ratios <- function(task, fold, learners, mtp, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)

  density_ratios <- matrix(nrow = nrow(natural$valid), ncol = task$time_horizon)
  fits <- vector("list", length = task$time_horizon)

  for (time in seq_len(task$time_horizon)) {
    i_train <- task$observed(natural$train, time - 1) %and% task$is_at_risk(natural$train, time)
    i_valid <- task$observed(natural$valid, time - 1) %and% task$is_at_risk(natural$valid, time)

    A_t <- current_trt(task$vars$A, time)

    # Treatment density ratio model — no longer includes censoring as a covariate
    trt_vars <- c("..i..lmtp_id", task$vars$history("A", time), A_t, "..i..lmtp_stack_indicator")
    stacked <- stack_data(natural$train, shifted$train, task$vars$A, task$vars$C, time)

    fit_trt <- run_ensemble(stacked[rep(i_train, 2), trt_vars], "..i..lmtp_stack_indicator",
                            learners, "binomial", "..i..lmtp_id",
                            control$.learners_trt_folds)

    # Separate censoring model
    if (!is.null(task$vars$C)) {
      cens_vars <- c("..i..lmtp_id", task$vars$history_cens(time), task$vars$C[time])
      fit_cens <- run_ensemble(natural$train[i_train, cens_vars], task$vars$C[time],
                               learners, "binomial", "..i..lmtp_id",
                               control$.learners_trt_folds)
      pred_cens <- rep(-999L, nrow(natural$valid))
      pred_cens[i_valid] <- predict(fit_cens, natural$valid[i_valid, ])
    } else {
      fit_cens <- NULL
      pred_cens <- 1
    }

    if (control$.return_full_fits) {
      fits[[time]] <- list(treatment = fit_trt, censoring = fit_cens)
    } else {
      fits[[time]] <- list(treatment = extract_sl_weights(fit_trt),
                           censoring = if (is.null(fit_cens)) NULL else extract_sl_weights(fit_cens))
    }

    pred <- rep(-999L, nrow(natural$valid))
    pred[i_valid] <- predict(fit_trt, natural$valid[i_valid, ])

    obs <- task$observed(natural$valid, time)
    at_risk <- task$is_at_risk(natural$valid, time)
    followed <- followed_rule(natural$valid, shifted$valid, A_t, mtp)

    pred <- ifelse(followed & !mtp, pmax(pred, 0.5), pred)
    density_ratios[, time] <- ((pred * obs * at_risk * followed) / (1 - pmin(pred, 0.999))) * (1 / pred_cens)

    pb()
  }

  list(ratios = density_ratios, fits = fits)
}

stack_data <- function(natural, shifted, trt, cens, time) {
  shifted_half <- natural

  if (length(trt) > 1 || time == 1) {
    shifted_half[, trt[[time]]] <- shifted[, trt[[time]]]
  }

  if (!is.null(cens)) {
    shifted_half[[cens[time]]] <- shifted[[cens[time]]]
  }

  out <- rbind(natural, shifted_half)
  out[["..i..lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(natural))
  out
}
