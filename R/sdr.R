cf_sdr <- function(Task, outcome, ratios, learners, lrnr_folds, pb) {
  out <- list()
  for (fold in seq_along(Task$folds)) {
    out[[fold]] <- future::future({
      estimate_sdr(
        get_folded_data(Task$natural, Task$folds, fold),
        get_folded_data(Task$shifted, Task$folds, fold),
        outcome, Task$node_list$outcome,
        Task$cens, Task$risk, Task$tau, Task$outcome_type,
        get_folded_data(ratios, Task$folds, fold)$train,
        learners, lrnr_folds, pb
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    natural = recombine_outcome(out, "natural", Task$folds),
    shifted = recombine_outcome(out, "shifted", Task$folds),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_sdr <- function(natural, shifted, outcome, node_list, cens, risk, tau,
                         outcome_type, ratios, learners, lrnr_folds, pb) {

  m_natural_train <- m_shifted_train <-
    cbind(matrix(nrow = nrow(natural$train), ncol = tau), natural$train[[outcome]])
  m_natural_valid <- m_shifted_valid <-
    cbind(matrix(nrow = nrow(natural$valid), ncol = tau), natural$valid[[outcome]])

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
      train_task <- sl3_task(natural$train[i & rt,], outcome, vars, outcome_type, "lmtp_id", lrnr_folds)

      natural_task <- sl3_task(natural$train[jt & rt,], NULL, vars, NULL, "lmtp_id")
      shift_task <- sl3_task(shifted$train[jt & rt,], NULL, vars, NULL, "lmtp_id")
      valid_natural_task <- sl3_task(natural$valid[jv & rv, ], NULL, vars, NULL, "lmtp_id")
      valid_shifted_task <- sl3_task(shifted$valid[jv & rv, ], NULL, vars, NULL, "lmtp_id")

      learners <- check_variation(natural$train[i & rt, ][[outcome]], learners)
      fit <- learners$train(train_task)
      fits[[t]] <- fit$fit_object
    }

    if (t < tau) {
      densratio <- transform_sdr(ratio_sdr(ratios, t, tau), t, tau, m_shifted_train, m_natural_train)
      natural$train[, pseudo] <- shifted$train[, pseudo] <- densratio

      train_task <- sl3_task(natural$train[i & rt,], pseudo, vars, "continuous", "lmtp_id", lrnr_folds)

      natural_task <- sl3_task(natural$train[jt & rt,], NULL, vars, NULL, "lmtp_id")
      shift_task <- sl3_task(shifted$train[jt & rt,], NULL, vars, NULL, "lmtp_id")
      valid_natural_task <- sl3_task(natural$valid[jv & rv, ], NULL, vars, NULL, "lmtp_id")
      valid_shifted_task <- sl3_task(shifted$valid[jv & rv, ], NULL, vars, NULL, "lmtp_id")

      learners <- check_variation(natural$train[i & rt, ][[pseudo]], learners)
      fit <- learners$train(train_task)
      fits[[t]] <- fit$fit_object
    }

    m_natural_train[jt & rt, t] <- bound(fit$predict(natural_task), 1e-05)
    m_shifted_train[jt & rt, t] <- bound(fit$predict(shift_task), 1e-05)
    m_natural_valid[jv & rv, t] <- bound(fit$predict(valid_natural_task), 1e-05)
    m_shifted_valid[jv & rv, t] <- bound(fit$predict(valid_shifted_task), 1e-05)

    m_natural_train[!rt, t] <- 0
    m_shifted_train[!rt, t] <- 0
    m_natural_valid[!rv, t] <- 0
    m_shifted_valid[!rv, t] <- 0

    pb()
  }

  list(
    natural = m_natural_valid,
    shifted = m_shifted_valid,
    fits = fits
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
