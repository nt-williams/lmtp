
# engine for density ratio estimation by classification
estimate_r <- function(data, A, C, shift, tau,
                       node_list, learner_stack = NULL, pb) {

  # setup
  n <- nrow(data)
  r <- list(natural = matrix(nrow = n, ncol = tau),
            shifted = matrix(nrow = n, ncol = tau))

  lapply(1:tau, function(t) {

    # progress bar
    progress_progress_bar(pb)

    # setup
    shifted <- shift_data(data, A[[t]], shift)
    d <- rbind(data, shifted)
    d$id <- rep(1:n, 2)
    d$shift_indicator <- c(rep(0, n), rep(1, n))
    task <- initiate_sl3_task(d, "shift_indicator", c(node_list[[t]]), "binomial", "id")
    ensemble <- initiate_ensemble("binomial", learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, task)

    # ratios
    pred <- bound(predict_sl3(fit, task), .Machine$double.eps)
    rat <- pred / (1 - truncate(pred))
    r$natural[, t] <<- rat[d$shift_indicator == 0] * C[, t]
    r$shifted[, t] <<- rat[d$shift_indicator == 1] * C[, t]
  })

  # returns
  return(r)
}

estimate_c <- function(data, C, Y, tau, node_list, learner_stack, cens) {

  # setup
  if (isFALSE(check_censoring(data, C, Y))) {
    lapply(1:tau, function(t) {
      cens[, t] <<- rep(1, nrow(cens))
    })
    return(cens)
  }

  no_cens <- data
  sapply(C, function(x) {
    no_cens[, x] <<- rep(1, nrow(no_cens))
  }, simplify = TRUE)

  lapply(1:tau, function(t) {

    browser()
    # setup
    if (t > 1) data <- data[data[[C[[t - 1]]]] == 1, ]
    fit_task <- initiate_sl3_task(data, C[[t]], node_list[[t]], "binomial", drop = TRUE)
    pred1_task <- initiate_sl3_task(no_cens, C[[t]], node_list[[t]], "binomial")
    # DO I also need the probabilities of being censored
    # then the weights for each observation... idk
    ensemble <- initiate_ensemble("binomial", learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, fit_task)

    # probability of not being censored
    pred <- bound(predict_sl3(fit, pred_task), .Machine$double.eps)
    rat <- pred / (1 - truncate(pred))
    cens[, t] <<- rat
  })

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
