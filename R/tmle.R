cf_tmle <- function(task, outcome, cumulated, ratios, learners, control, pb) {
  out <- vector("list", length = length(task$folds))

  if (isFALSE(cumulated)) {
    ratios <- matrix(t(apply(ratios, 1, cumprod)),
                     nrow = nrow(ratios),
                     ncol = ncol(ratios))
  }

  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_tmle(
        get_folded_data(task$natural, task$folds, fold),
        get_folded_data(task$shifted[, unlist(task$trt), drop = F], task$folds, fold),
        task$trt,
        outcome,
        task$node_list$outcome,
        task$cens,
        task$risk,
        task$tau,
        task$conditional[task$folds[[fold]]$training_set, , drop = FALSE],
        task$outcome_type,
        get_folded_data(ratios, task$folds, fold)$train,
        task$weights[task$folds[[fold]]$training_set],
        learners,
        control,
        pb
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    natural = recombine_outcome(out, "natural", task$folds),
    shifted = cbind(recombine_outcome(out, "shifted", task$folds), task$natural[["tmp_lmtp_scaled_outcome"]]),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_tmle <- function(natural, shifted, trt, outcome, node_list, cens,
                          risk, tau, conditional, outcome_type, ratios, weights, learners, control, pb) {
  m_natural_train <- m_shifted_train <- matrix(nrow = nrow(natural$train), ncol = tau)
  m_natural_valid <- m_shifted_valid <- matrix(nrow = nrow(natural$valid), ncol = tau)

  fits <- vector("list", length = tau)
  for (t in tau:1) {
    i  <- censored(natural$train, cens, t)$i
    jt <- censored(natural$train, cens, t)$j
    jv <- censored(natural$valid, cens, t)$j
    rt <- at_risk(natural$train, risk, t)
    rv <- at_risk(natural$valid, risk, t)

    in_conditioning_set <- as.logical(apply(conditional[, (t + 1):(tau + 1), drop = FALSE], 1, prod))

    pseudo <- paste0("tmp_lmtp_pseudo", t)
    vars <- node_list[[t]]

    if (t != tau) {
      outcome <- paste0("tmp_lmtp_pseudo", t + 1)
      outcome_type <- "continuous"
    }

    fit <- run_ensemble(natural$train[i & rt & in_conditioning_set, c("lmtp_id", vars, outcome)],
                        outcome,
                        learners,
                        outcome_type,
                        "lmtp_id",
                        control$.learners_outcome_folds)

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

    under_shift_train <- natural$train[jt & rt, c("lmtp_id", vars)]
    under_shift_train[, trt_t] <- shifted$train[jt & rt, trt_t]

    under_shift_valid <- natural$valid[jv & rv, c("lmtp_id", vars)]
    under_shift_valid[, trt_t] <- shifted$valid[jv & rv, trt_t]

    m_natural_train[jt & rt, t] <- bound(SL_predict(fit, natural$train[jt & rt, c("lmtp_id", vars)]), 1e-05)
    m_shifted_train[jt & rt, t] <- bound(SL_predict(fit, under_shift_train), 1e-05)
    m_natural_valid[jv & rv, t] <- bound(SL_predict(fit, natural$valid[jv & rv, c("lmtp_id", vars)]), 1e-05)
    m_shifted_valid[jv & rv, t] <- bound(SL_predict(fit, under_shift_valid), 1e-05)

    wts <- (ratios[, t] * weights)[i & rt & in_conditioning_set]

    fit <- sw(
      glm(
        natural$train[i & rt & in_conditioning_set, ][[outcome]] ~ offset(qlogis(m_natural_train[i & rt & in_conditioning_set, t])),
        weights = wts,
        family = "binomial"
      )
    )

    natural$train[jt & rt, pseudo] <- bound(plogis(qlogis(m_shifted_train[jt & rt, t]) + coef(fit)))
    m_natural_valid[jv & rv, t] <- bound(plogis(qlogis(m_natural_valid[jv & rv, t]) + coef(fit)))
    m_shifted_valid[jv & rv, t] <- bound(plogis(qlogis(m_shifted_valid[jv & rv, t]) + coef(fit)))

    natural$train[!rt, pseudo] <- 0
    m_natural_valid[!rv, t] <- 0
    m_shifted_valid[!rv, t] <- 0

    pb()
  }

  list(natural = m_natural_valid,
       shifted = m_shifted_valid,
       fits = fits)
}
