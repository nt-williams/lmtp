
# engine for density ratio estimation by classification
estimate_r <- function(data, A, C, shift, tau,
                       node_list, learner_stack = NULL, pb) {

  # global setup
  n <- nrow(data)
  r <- list(natural = matrix(nrow = n, ncol = tau),
            shifted = matrix(nrow = n, ncol = tau))

  lapply(1:tau, function(t) {

    # progress bar
    progress_progress_bar(pb)

    # setup
    i <- !is.na(data[[A[[t]]]])
    shifted <- shift_data(data, A[[t]], shift)
    d <- rbind(data, shifted)
    d$id <- rep(1:n, 2)
    d$shift_indicator <- c(rep(0, n), rep(1, n))
    d <- subset(d, rep(i, 2))
    task <- initiate_sl3_task(d, "shift_indicator", node_list[[t]], "binomial", "id")
    ensemble <- initiate_ensemble("binomial", learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, task)

    # ratios
    pred <- bound(predict_sl3(fit, task), .Machine$double.eps)
    rat <- pred / (1 - truncate(pred))
    r$natural[!i, t] <<- 0
    r$shifted[!i, t] <<- 0
    r$natural[i, t] <<- rat[d$shift_indicator == 0] * C[i, t]
    r$shifted[i, t] <<- rat[d$shift_indicator == 1] * C[i, t]
  })

  # returns
  return(r)
}

# engine for estimation of censoring mechanism
estimate_c <- function(data, C, Y, tau, node_list, learner_stack) {

  # global setup
  cens <- check_censoring(data, C, Y, tau)

  if (all(is.na(cens))) {
    lapply(1:tau, function(t) {

      # setup
      fit_task <- suppressWarnings(initiate_sl3_task(data, C[[t]], node_list[[t]], "binomial", drop = TRUE))
      ensemble <- initiate_ensemble("binomial", learner_stack)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # probability of not being censored
      pred <- bound(predict_sl3(fit, fit_task), .Machine$double.eps)
      rat <- mean(data[, C[[t]]]) / pred
      cens[, t] <<- rat
    })
  }

  # returns
  return(cens)
}

use_dens_ratio <- function(r, tau, n, max, what) {
  switch(
    what,
    "tmle" = ratio_ite(r = r, tau = tau, n = n),
    "ipw" = ratio_ite(r = r, tau = tau, n = n),
    "eif" = ratio_ite(r = r, tau = tau, n = n),
    "sdr" = ratio_sdr(r = r, tau = tau, max = max)
  )
}

ratio_ite <- function(r, tau, n) {
  out <- matrix(t(apply(r$natural, 1, cumprod)), nrow = n, ncol = tau)
  return(out)
}

ratio_sdr <- function(r, tau, max) {
  out <- t(apply(r$natural[, (tau + 1):max, drop = FALSE], 1, cumprod))
  if (tau == max - 1) out <- t(out)
  return(out)
}
