cf_sub <- function(task, outcome, learners, control, pb) {
  out <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_sub(
        get_folded_data(task$natural, task$folds, fold),
        get_folded_data(task$shifted[, unlist(task$trt), drop = F], task$folds, fold),
        task$trt,
        outcome,
        task$node_list$outcome,
        task$cens,
        task$risk,
        task$tau,
        task$outcome_type,
        learners,
        control,
        pb
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    m = recombine_outcome(out, "m", task$folds),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_sub <- function(natural, shifted, trt, outcome, node_list, cens, risk,
                         tau, outcome_type, learners, control, pb) {

  m <- matrix(nrow = nrow(natural$valid), ncol = tau)
  fits <- vector("list", length = tau)

  for (t in tau:1) {
    i  <- censored(natural$train, cens, t)$i
    jt <- censored(natural$train, cens, t)$j
    jv <- censored(natural$valid, cens, t)$j
    rt <- at_risk(natural$train, risk, t)
    rv <- at_risk(natural$valid, risk, t)

    pseudo <- paste0("tmp_lmtp_pseudo", t)
    vars <- node_list[[t]]

    if (t != tau) {
      outcome <- paste0("tmp_lmtp_pseudo", t + 1)
      outcome_type <- "continuous"
    }

    learners <- check_variation(natural$train[i & rt, ][[outcome]], learners)

    fit <- run_ensemble(
      natural$train[i & rt, ][[outcome]],
      natural$train[i & rt, vars],
      learners,
      outcome_type,
      id = natural$train[i & rt, ][["lmtp_id"]],
      control$.learners_outcome_folds
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    if (length(trt) > 1) {
      trt_t <- trt[[t]]
    } else {
      trt_t <- trt[[1]]
    }

    under_shift_train <- natural$train[jt & rt, vars]
    under_shift_train[, trt_t] <- shifted$train[jt & rt, trt_t]

    under_shift_valid <- natural$valid[jv & rv, vars]
    under_shift_valid[, trt_t] <- shifted$valid[jv & rv, trt_t]

    natural$train[jt & rt, pseudo] <- bound(SL_predict(fit, under_shift_train), 1e-05)
    m[jv & rv, t] <- bound(SL_predict(fit, under_shift_valid), 1e-05)

    natural$train[!rt, pseudo] <- 0
    m[!rv, t] <- 0

    pb()
  }

  list(m = m, fits = fits)
}
