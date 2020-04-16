
# engine for density ratio estimation by classification
estimate_r <- function(data, trt, cens, C, shift, tau,
                       node_list, learners = NULL, pb) {

  # global setup
  n <- nrow(data)
  r <- list(natural = matrix(nrow = n, ncol = tau),
            shifted = matrix(nrow = n, ncol = tau))

  if (!is.null(shift)) {

    for (t in 1:tau) {
      # progress bar
      progress_progress_bar(pb)

      # setup
      i        <- create_censoring_indicators(data, cens, tau)$i
      shifted  <- shift_data(data, trt[[t]], shift)
      d        <- rbind(data, shifted)
      d$id     <- rep(1:n, 2)
      d$si     <- c(rep(0, n), rep(1, n))
      d        <- subset(d, rep(i, 2))
      task     <- initiate_sl3_task(d, "si", node_list[[t]], "binomial", "id")
      ensemble <- initiate_ensemble("binomial", learners)

      # run SL
      fit <- run_ensemble(ensemble, task)

      # ratios
      pred             <- bound(predict_sl3(fit, task), .Machine$double.eps)
      rat              <- pred / (1 - truncate(pred))
      r$natural[!i, t] <- 0
      r$shifted[!i, t] <- 0
      r$natural[i, t]  <- rat[d$si == 0] * C[i, t]
      r$shifted[i, t]  <- rat[d$si == 1] * C[i, t]
    }

  } else {

    for (t in 1:tau) {
      # progress bar
      progress_progress_bar(pb)

      # setup
      i <- create_censoring_indicators(data, cens, tau)$i

      # propensity
      r$natural[!i, t] <- 0
      r$shifted[!i, t] <- 0
      r$natural[i, t]  <- C[i, t]
      r$shifted[i, t]  <- C[i, t]
    }

  }

  # returns
  return(r)
}

# engine for estimation of censoring mechanism
estimate_c <- function(data, C, outcome, tau, node_list, learners) {

  # global setup
  cens <- check_censoring(data, C, outcome, tau)

  if (all(is.na(cens))) {
    for (t in 1:tau) {
      # setup
      fit_task <- suppressWarnings(initiate_sl3_task(data, C[[t]], node_list[[t]], "binomial", drop = TRUE))
      ensemble <- initiate_ensemble("binomial", learners)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # probability of not being censored
      cens[, t] <- mean(data[, C[[t]]]) /
        bound(predict_sl3(fit, fit_task), .Machine$double.eps)
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
