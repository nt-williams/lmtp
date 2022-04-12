cf_sub <- function(Task, outcome, learners, lrnr_folds, pb) {
  out <- list()

  for (fold in seq_along(Task$folds)) {
    out[[fold]] <- future::future({
      estimate_sub(
        get_folded_data(Task$natural, Task$folds, fold),
        get_folded_data(Task$shifted, Task$folds, fold),
        outcome,
        Task$node_list$outcome, Task$cens,
        Task$risk, Task$tau, Task$outcome_type,
        learners, lrnr_folds, pb
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    m = recombine_outcome(out, "m", Task$folds),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_sub <- function(natural, shifted, outcome, node_list, cens, risk,
                         tau, outcome_type, learners, lrnr_folds, pb) {

  m <- matrix(nrow = nrow(natural$valid), ncol = tau)
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

    train_task <- sl3_task(natural$train[i & rt,], outcome, vars, outcome_type, "lmtp_id", lrnr_folds)

    shift_task <- sl3_task(shifted$train[jt & rt,], NULL, vars, NULL, "lmtp_id")
    valid_shift_task <- sl3_task(shifted$valid[jv & rv, ], NULL, vars, NULL, "lmtp_id")

    learners <- check_variation(natural$train[i & rt, ][[outcome]], learners)
    fit <- learners$train(train_task)
    fits[[t]] <- fit$fit_object

    natural$train[jt & rt, pseudo] <- bound(fit$predict(shift_task), 1e-05)
    m[jv & rv, t] <- bound(fit$predict(valid_shift_task), 1e-05)

    natural$train[!rt, pseudo] <- 0
    m[!rv, t] <- 0

    pb()
  }

  list(m = m, fits = fits)
}
