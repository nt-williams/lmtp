
# the engine for density ratio nuisance parameter
# currently only works for a single time point
# estimate_r_sl <- function(data, A, shifted_down, tau = NULL,
#                           node_list = NULL, r = NULL, g = NULL,
#                           learner_stack = NULL) {
#
#   # setup
#   ensemble <- initiate_ensemble(data, A, c("W1", "W2"), "density", learner_stack)
#   to_predict <- initiate_sl3_task(shifted_down, A, c("W1", "W2"), "density")
#
#   # run SL
#   fit <- run_ensemble(ensemble)
#
#   # predict on shifted data
#   g[, 1] <- predict_sl3_density(fit, ensemble$task)
#   g[, 2] <- predict_sl3_density(fit, to_predict)
#
#   # density ratios
#   r[, 1] <- g[, 2] / g[, 1]
#
#   # returns
#   return(r)
# }

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
    pred <- bound(predict_sl3_nondensity(fit, task), .Machine$double.eps)
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


