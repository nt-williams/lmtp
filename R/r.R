
# the engine for density ratio nuisance parameter
# currently only works for a single time point
estimate_r_sl <- function(data, A, shifted_down, tau = NULL,
                          node_list = NULL, r = NULL, g = NULL,
                          learner_stack = NULL) {

  # setup
  ensemble <- initiate_ensemble(data, A, c("W1", "W2"), "density", learner_stack)
  to_predict <- initiate_sl3_task(shifted_down, A, c("W1", "W2"), "density")

  # run SL
  fit <- run_ensemble(ensemble)

  # predict on shifted data
  g[, 1] <- predict_sl3_density(fit, ensemble$task)
  g[, 2] <- predict_sl3_density(fit, to_predict)

  # density ratios
  r[, 1] <- g[, 2] / g[, 1]

  # returns
  return(r)
}
