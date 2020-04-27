
# engine for density ratio estimation by classification
estimate_r <- function(training, validation, trt, cens, C,
                       shift, tau, node_list, learners = NULL, pb) {

  # global setup
  nt <- nrow(training)
  nv <- nrow(validation)
  r  <- list(natural = matrix(nrow = nv, ncol = tau),
             shifted = matrix(nrow = nv, ncol = tau))

  if (!is.null(shift)) {

    for (t in 1:tau) {
      # progress bar
      progress_progress_bar(pb)

      # setup
      train_stck <- prepare_r_engine(training, shift_data(training, trt[[t]], shift), nt)
      valid_stck <- prepare_r_engine(validation, shift_data(validation, trt[[t]], shift), nt)

      # create sl3 tasks for training and validation sets
      fit_task  <-
        initiate_sl3_task(subset(train_stck, rep(
          create_censoring_indicators(training, cens, tau)$i, 2
        )), "si", node_list[[t]], "binomial", "id")
      pred_task <- suppressWarnings(initiate_sl3_task(valid_stck, "si", node_list[[t]], "binomial", "id"))
      ensemble  <- initiate_ensemble("binomial", learners)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # ratios
      pred           <- bound(predict_sl3(fit, pred_task), .Machine$double.eps)
      rat            <- pred / (1 - truncate(pred))
      r$natural[, t] <- rat[valid_stck$si == 0] * C[, t]
      r$shifted[, t] <- rat[valid_stck$si == 1] * C[, t]
    }

  } else {

    for (t in 1:tau) {
      # progress bar
      progress_progress_bar(pb)

      # propensity
      r$natural[, t] <- C[, t]
      r$shifted[, t] <- C[, t]
    }

  }

  r$natural <- check_extreme_ratio(r$natural)
  r$shifted <- check_extreme_ratio(r$shifted)

  # returns
  return(r)
}

# engine for estimation of censoring mechanism
estimate_c <- function(data, training, validation, C,
                       outcome, tau, node_list, learners) {

  # global setup
  cens <- check_censoring(data, validation, C, outcome, tau)

  if (all(is.na(cens))) {
    for (t in 1:tau) {
      # setup
      fit_task  <- suppressWarnings(initiate_sl3_task(training, C[[t]], node_list[[t]], "binomial", drop = TRUE))
      pred_task <- suppressWarnings(initiate_sl3_task(validation, C[[t]], node_list[[t]], "binomial", drop = TRUE))
      ensemble  <- initiate_ensemble("binomial", learners)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # probability of not being censored
      cens[, t] <- mean(data[, C[[t]]]) /
        bound(predict_sl3(fit, pred_task), .Machine$double.eps)
    }
  }

  # returns
  return(cens)
}

use_dens_ratio <- function(ratio, tau, n, max_tau, what_estim) {
  switch(
    what_estim,
    "tml" = ratio_ite(ratio = ratio, tau = tau, n = n),
    "ipw" = ratio_ite(ratio = ratio, tau = tau, n = n),
    "eif" = ratio_ite(ratio = ratio, tau = tau, n = n),
    "sdr" = ratio_sdr(ratio = ratio, tau = tau, max_tau = max_tau)
  )
}

ratio_ite <- function(ratio, tau, n) {
  out <- matrix(t(apply(ratio$natural, 1, cumprod)), nrow = n, ncol = tau)
  return(out)
}

ratio_sdr <- function(ratio, tau, max_tau) {
  out <- t(apply(ratio$natural[, (tau + 1):max_tau, drop = FALSE], 1, cumprod))
  if (tau == max_tau - 1) out <- t(out)
  return(out)
}
