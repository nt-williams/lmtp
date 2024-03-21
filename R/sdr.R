cf_sdr <- function(task, outcome, ratios, learners, control, progress_bar) {
  out <- list()
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_sdr(get_folded_data(task$natural, task$folds, fold),
                   get_folded_data(task$shifted[, task$trt, drop = F], task$folds, fold),
                   outcome, task$node_list$outcome,
                   task$cens,
                   task$risk,
                   task$id,
                   task$tau,
                   task$outcome_type,
                   get_folded_data(ratios, task$folds, fold)$train,
                   learners,
                   control,
                   progress_bar)
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(natural = recombine_outcome(out, "natural", task$folds),
       shifted = recombine_outcome(out, "shifted", task$folds),
       fits = lapply(out, function(x) x[["fits"]]))
}

estimate_sdr <- function(natural, shifted, outcome, node_list, cens, risk, id, tau,
                         outcome_type, ratios, learners, control, progress_bar) {

  m_natural_train <- m_shifted_train <- cbind(matrix(nrow = nrow(natural$train),
                                                     ncol = tau),
                                              natural$train[[outcome]])

  m_natural_valid <- m_shifted_valid <- cbind(matrix(nrow = nrow(natural$valid),
                                                     ncol = tau),
                                              natural$valid[[outcome]])

  fits <- list()

  for (t in tau:1) {
    i  <- censored(natural$train, cens, t)$i
    jt <- censored(natural$train, cens, t)$j
    jv <- censored(natural$valid, cens, t)$j
    rt <- at_risk(natural$train, risk, t)
    rv <- at_risk(natural$valid, risk, t)

    pseudo <- paste0("tmp_lmtp_pseudo", t)
    vars <- node_list[[t]]

    if (t == tau) {
      learners <- check_variation(natural$train[i & rt, ][[outcome]], learners)

      fit <- run_ensemble(natural$train[i & rt, c(id, vars, outcome)],
                          outcome,
                          learners,
                          outcome_type,
                          id,
                          control$.learners_outcome_folds)
    }

    if (t < tau) {
      densratio <- transform_sdr(compute_weights(ratios, t + 1, tau),
                                 t,
                                 tau,
                                 m_shifted_train,
                                 m_natural_train)

      natural$train[, pseudo] <- shifted$train[, pseudo] <- densratio

      learners <- check_variation(natural$train[i & rt, ][[pseudo]], learners)

      fit <- run_ensemble(natural$train[i & rt, c(id, vars, pseudo)],
                          pseudo,
                          learners,
                          "continuous",
                          id,
                          control$.learners_outcome_folds)
    }

<<<<<<< HEAD
    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

=======
>>>>>>> 3d968bb (Fix for issue #130)
    trt_var <- names(shifted$train)[t]
    under_shift_train <- natural$train[jt & rt, vars]
    under_shift_train[[trt_var]] <- shifted$train[jt & rt, trt_var]

    under_shift_valid <- natural$valid[jv & rv, vars]
    under_shift_valid[[trt_var]] <- shifted$valid[jv & rv, trt_var]

    m_natural_train[jt & rt, t] <- bound(SL_predict(fit, natural$train[jt & rt, vars]), 1e-05)
    m_shifted_train[jt & rt, t] <- bound(SL_predict(fit, under_shift_train), 1e-05)
    m_natural_valid[jv & rv, t] <- bound(SL_predict(fit, natural$valid[jv & rv, vars]), 1e-05)
    m_shifted_valid[jv & rv, t] <- bound(SL_predict(fit, under_shift_valid), 1e-05)

    m_natural_train[!rt, t] <- 0
    m_shifted_train[!rt, t] <- 0
    m_natural_valid[!rv, t] <- 0
    m_shifted_valid[!rv, t] <- 0

    progress_bar()
  }

  list(natural = m_natural_valid,
       shifted = m_shifted_valid,
       fits = fits)
}
