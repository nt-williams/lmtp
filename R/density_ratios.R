cf_r <- function(task, learners, mtp, control, pb) {
  ans <- vector("list", length = length(task$folds))

  if (length(learners) == 1 && learners == "mean") {
    warning("Using 'mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('glm', 'cv_glmnet')`.",
            call. = FALSE)
  }

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_r(task, fold, learners, mtp, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  ans <- list(ratios = recombine(rbind_depth(ans, "ratios"), task$folds),
              fits = lapply(ans, function(x) x[["fits"]]))

  ans$ratios <- trim(ans$ratios, control$.trim)
  ans
}

estimate_r <- function(task, fold, learners, mtp, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)

  density_ratios <- matrix(nrow = nrow(natural$valid), ncol = task$tau)
  fits <- vector("list", length = task$tau)

  for (t in 1:task$tau) {
    i <- ii(task$observed(natural$train, t - 1), task$R(natural$train, t))
    i <- rep(i, 2)

    if (length(task$vars$A) > 1) {
      A_t <- task$vars$A[[t]]
    } else {
      A_t <- task$vars$A[[1]]
    }

    vars <- c("..i..lmtp_id", task$vars$history("A", t), A_t, task$vars$C[t], "..i..lmtp_stack_indicator")
    stacked <- stack_data(natural$train, shifted$train, task$vars$A, task$vars$C, t)

    fit <- run_ensemble(stacked[i, vars], "..i..lmtp_stack_indicator",
                        learners, "binomial", "..i..lmtp_id",
                        control$.learners_trt_folds)

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    i <- ii(task$observed(natural$valid, t - 1), task$R(natural$valid, t))

    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    pred[i, ] <- predict(fit, natural$valid[i, ])

    obs <- task$observed(natural$valid, t)
    at_risk <- task$R(natural$valid, t)
    followed <- followed_rule(natural$valid, shifted$valid, A_t, mtp)

    pred <- ifelse(followed & !mtp, pmax(pred, 0.5), pred)
    density_ratios[, t] <- (pred * obs * at_risk * followed) / (1 - pmin(pred, 0.999))

    pb()
  }

  list(ratios = density_ratios, fits = fits)
}

stack_data <- function(natural, shifted, A, C, t) {
  shifted_half <- natural

  if (length(A) > 1 || t == 1) {
    shifted_half[, A[[t]]] <- shifted[, A[[t]]]
  }

  if (!is.null(C)) {
    shifted_half[[C[t]]] <- shifted[[C[t]]]
  }

  out <- rbind(natural, shifted_half)
  out[["..i..lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(natural))
  out
}
