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
    i <- task$observed(natural$train, time - 1) %and% task$is_at_risk(natural$train, time)
    i <- rep(i, 2)

    A_t <- current_trt(task$vars$A, time)

    vars <- c("..i..lmtp_id", task$vars$history("A", time), A_t, task$vars$C[time], "..i..lmtp_stack_indicator")
    stacked <- stack_data(natural$train, shifted$train, task$vars$A, task$vars$C, time)

    fit <- run_ensemble(stacked[i, vars], "..i..lmtp_stack_indicator",
                        learners, "binomial", "..i..lmtp_id",
                        control$.learners_trt_folds)

    if (control$.return_full_fits) {
      fits[[time]] <- fit
    } else {
      fits[[time]] <- extract_sl_weights(fit)
    }

    i <- task$observed(natural$valid, time - 1) %and% task$is_at_risk(natural$valid, time)

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred[i, ] <- predict(fit, natural$valid[i, ])

    obs <- task$observed(natural$valid, time)
    at_risk <- task$is_at_risk(natural$valid, time)
    followed <- followed_rule(natural$valid, shifted$valid, A_t, mtp)

    pred <- ifelse(followed & !mtp, pmax(pred, 0.5), pred)
    density_ratios[, time] <- (pred * obs * at_risk * followed) / (1 - pmin(pred, 0.999))

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
