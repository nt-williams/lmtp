cf_tmle <- function(Task, outcome, ratios, learners, lrnr_folds, full_fits, pb) {
  out <- list()

  ratios <- matrix(
    t(apply(ratios, 1, cumprod)),
    nrow = nrow(ratios),
    ncol = ncol(ratios)
  )

  for (fold in seq_along(Task$folds)) {
    out[[fold]] <- future::future({
      estimate_tmle(
        get_folded_data(Task$natural, Task$folds, fold),
        get_folded_data(Task$shifted, Task$folds, fold),
        outcome, Task$node_list$outcome, Task$cens, Task$risk,
        Task$tau, Task$outcome_type,
        get_folded_data(ratios, Task$folds, fold)$train,
        Task$weights[Task$folds[[fold]]$training_set],
        learners, lrnr_folds, pb, full_fits
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    natural = recombine_outcome(out, "natural", Task$folds),
    shifted = cbind(recombine_outcome(out, "shifted", Task$folds), Task$natural[["tmp_lmtp_scaled_outcome"]]),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_tmle <- function(natural, shifted, outcome, node_list, cens, risk, tau, outcome_type, ratios, weights, learners, lrnr_folds, pb, full_fits) {
  m_natural_train <- m_shifted_train <- matrix(nrow = nrow(natural$train), ncol = tau)
  m_natural_valid <- m_shifted_valid <- matrix(nrow = nrow(natural$valid), ncol = tau)

  fits <- list()
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
      id = natural$train[i & rt,][["lmtp_id"]],
      lrnr_folds
    )

    if (full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    m_natural_train[jt & rt, t] <- bound(SL_predict(fit, natural$train[jt & rt, vars]), 1e-05)
    m_shifted_train[jt & rt, t] <- bound(SL_predict(fit, shifted$train[jt & rt, vars]), 1e-05)
    m_natural_valid[jv & rv, t] <- bound(SL_predict(fit, natural$valid[jv & rv, vars]), 1e-05)
    m_shifted_valid[jv & rv, t] <- bound(SL_predict(fit, shifted$valid[jv & rv, vars]), 1e-05)

    wts <- {
      if (is.null(weights))
        ratios[i & rt, t]
      else
        ratios[i & rt, t] * weights[i & rt]
    }

    fit <- sw(
      glm(
        natural$train[i & rt, ][[outcome]] ~ offset(qlogis(m_natural_train[i & rt, t])),
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

  list(
    natural = m_natural_valid,
    shifted = m_shifted_valid,
    fits = fits
  )
}
