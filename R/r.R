
# engine for density ratio estimation by classification
estimate_r_sl <- function(data, A, shift, tau,
                          node_list, learner_stack = NULL) {

  # setup
  n <- nrow(data)
  r <- list(natural = matrix(nrow = n, ncol = tau),
            shifted = matrix(nrow = n, ncol = tau))

  lapply(1:tau, function(t) {
    # setup
    shifted <- shift_data(data, A[[t]], shift)
    d <- rbind(data, shifted)
    d$id <- rep(1:n, 2)
    d$shift_indicator <- c(rep(0, n), rep(1, n))
    task <- initiate_sl3_task(d, "shift_indicator", node_list[[t]], "binomial", "id")
    ensemble <- initiate_ensemble("binomial", learner_stack)

    # run SL
    fit <- run_ensemble(ensemble, task)

    # ratios
    pred <- bound(predict_sl3(fit, task), .Machine$double.eps)
    rat <- pred / (1 - truncate(pred))
    r$natural[, t] <<- rat[d$shift_indicator == 0]
    r$shifted[, t] <<- rat[d$shift_indicator == 1]
  })

  rn <- matrix(t(apply(r$natural, 1, cumprod)), nrow = n, ncol = tau)
  rd <- rn / (r$natural * r$shifted)

  # returns
  out <- list(rn = rn,
              rd = rd)

  return(out)
}


