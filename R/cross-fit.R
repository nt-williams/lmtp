
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
    out[[i]] <- future::future({
      estimate_c(data, folded[[i]]$train, folded[[i]]$valid,
                 C, outcome, tau, node_list, learners)
      }, packages = "lmtp")
  }
  out <- get_future_result(out)
  return(out)
}

cf_r <- function(data, shift, V, trt, cens, C,
                 tau, node_list, learners, pb) {

  out <- list()
  for (i in 1:V) {
    out[[i]] <- future::future({
      estimate_r(data[[i]]$train, data[[i]]$valid, trt,
                           cens, C[[i]], shift, tau, node_list, learners, pb)
    }, packages = "lmtp")
  }
  out <- get_future_result(out)
  return(out)
}

cf_sub <- function(data, shifted, V, outcome, node_list, C, tau,
                   outcome_type, learners, m, pb) {

  out <- list()
  for (i in 1:V) {
    out[[i]] <- future::future({
      estimate_sub(data[[i]]$train, shifted[[i]]$train, shifted[[i]]$valid,
                   outcome, node_list, C, tau, outcome_type,
                   learners, m[[i]]$valid, pb)
    }, packages = "lmtp")
  }
  out <- get_future_result(out)
  return(Reduce(rbind, out))

}

cf_tmle <- function(data, shifted, V, outcome, node_list, C, tau,
                    outcome_type, m_natural, m_shifted, r, learners, pb) {

  m <- list()
  for (i in 1:V) {
    m[[i]] <- future::future({
      estimate_tmle(data[[i]]$train, shifted[[i]]$train, data[[i]]$valid,
                    shifted[[i]]$valid, outcome, node_list, C, tau,
                    outcome_type, m_natural[[i]], m_shifted[[i]], r[[i]],
                    learners, pb)
    }, packages = "lmtp")
  }
  m <- get_future_result(m)
  out <- list(natural = Reduce(rbind, lapply(m, function(x) x[["natural"]])),
              shifted = Reduce(rbind, lapply(m, function(x) x[["shifted"]])))
  return(out)
}

cf_sdr <- function(data, shifted, V, outcome, node_list, C, tau,
                   outcome_type, m_natural, m_shifted, r, learners, pb) {
  m <- list()
  for (i in 1:V) {
    m[[i]] <- future::future({
      estimate_sdr(data[[i]]$train, shifted[[i]]$train,data[[i]]$valid,
                   shifted[[i]]$valid,outcome, node_list, C, tau, tau,
                   outcome_type, learners,m_natural[[i]], m_shifted[[i]],
                   r[[i]], pb)
    }, packages = "lmtp")
  }
  m <- get_future_result(m)
  out <- list(natural = Reduce(rbind, lapply(m, function(x) x[["natural"]])),
              shifted = Reduce(rbind, lapply(m, function(x) x[["shifted"]])))
  return(out)
}
