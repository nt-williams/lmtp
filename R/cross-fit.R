
setup_cv <- function(data, V = 10) {
  out <- origami::make_folds(data, V = V)
  return(out)
}

get_folded_data <- function(data, folds) {
  out <- list()
  for (i in 1:length(folds)) {
    out[[i]] <- list()
    out[[i]][["train"]] <- data[folds[[i]]$training_set, ]
    out[[i]][["valid"]] <- data[folds[[i]]$validation_set, ]
  }
  return(out)
}

cf_cens <- function(data, folded, V, C, outcome, tau, node_list, learners) {
  out <- list()
  for (i in 1:V) {
    out[[i]] <- estimate_c(data, folded[[i]]$train, folded[[i]]$valid, C,
                           outcome, tau, node_list, learners)
  }
  return(out)
}

cf_r <- function(data, shift, V, trt, cens, C,
                 tau, node_list, learners, pb) {

  out <- list()
  for (i in 1:V) {
    out[[i]] <- estimate_r(data[[i]]$train, data[[i]]$valid, trt,
                           cens, C[[i]], shift, tau, node_list, learners, pb)
  }
  return(out)
}

recombine_ipw <- function(r) {
  return(Reduce(rbind, Reduce(rbind, r)[, "natural"]))
}

cf_sub <- function(data, shifted, V, outcome, node_list, C, tau,
                   outcome_type, learners, m, pb) {

  out <- list()
  for (i in 1:V) {
    out[[i]] <-
      estimate_sub(data[[i]]$train, shifted[[i]]$train, shifted[[i]]$valid,
                   outcome, node_list, C, tau, outcome_type,
                   learners, m[[i]]$valid, pb)
  }
  return(Reduce(rbind, out))

}
