
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

cf_cens <- function(data, folded, V, C, outcome, tau,
                    node_list, learners, weights_c) {
  out <- list()
  fopts <- options("lmtp.bound")
  for (i in 1:V) {
    out[[i]] <- future::future({
      options(fopts)
      estimate_c(data, folded[[i]]$train, folded[[i]]$valid,
                 C, outcome, tau, node_list, learners, weights_c[[i]])
      }, packages = "lmtp")
  }
  out <- future::values(out)
  return(out)
}

cf_r <- function(data, shift, V, trt, cens, C, tau,
                 node_list, learners, pb, weights_r) {
  fopts <- options("lmtp.bound")
  out <- list()
  for (i in 1:V) {
    out[[i]] <- future::future({
      options(fopts)
      estimate_r(data[[i]]$train, data[[i]]$valid, trt, cens, C[[i]],
                 shift, tau, node_list, learners, pb, weights_r[[i]])
    }, packages = "lmtp")
  }
  out <- future::values(out)
  return(out)
}

cf_sub <- function(data, shifted, V, outcome, node_list, C, tau,
                   outcome_type, learners, m, pb) {
  fopts <- options("lmtp.bound")
  out <- list()
  for (i in 1:V) {
    out[[i]] <- future::future({
      options(fopts)
      estimate_sub(data[[i]]$train, shifted[[i]]$train, shifted[[i]]$valid,
                   outcome, node_list, C, tau, outcome_type,
                   learners, m[[i]]$valid, pb)
    }, packages = "lmtp")
  }
  out <- future::values(out)
  return(Reduce(rbind, out))

}

cf_tmle <- function(data, shifted, V, outcome, node_list, C, tau, outcome_type,
                    m_natural, m_shifted, r, learners, pb, weights_m) {

  fopts <- options("lmtp.bound")
  m <- list()
  for (i in 1:V) {
    m[[i]] <- future::future({
      options(fopts)
      estimate_tmle(data[[i]]$train, shifted[[i]]$train, data[[i]]$valid,
                    shifted[[i]]$valid, outcome, node_list, C, tau,
                    outcome_type, m_natural[[i]], m_shifted[[i]], r[[i]],
                    learners, pb, weights_m[[i]])
    }, packages = "lmtp")
  }
  m <- future::values(m)
  out <- list(natural = Reduce(rbind, lapply(m, function(x) x[["natural"]])),
              shifted = Reduce(rbind, lapply(m, function(x) x[["shifted"]])),
              sl_weights = lapply(m, function(x) x[["sl_weights"]]))
  return(out)
}

cf_sdr <- function(data, shifted, V, outcome, node_list, C, tau, outcome_type,
                   m_natural, m_shifted, r, learners, pb, weights_m) {

  fopts <- options("lmtp.bound")
  m <- list()
  for (i in 1:V) {
    m[[i]] <- future::future({
      options(fopts)
      estimate_sdr(data[[i]]$train, shifted[[i]]$train,data[[i]]$valid,
                   shifted[[i]]$valid,outcome, node_list, C, tau, tau,
                   outcome_type, learners,m_natural[[i]], m_shifted[[i]],
                   r[[i]], pb, weights_m[[i]])
    }, packages = "lmtp")
  }
  m <- future::values(m)
  out <- list(natural = Reduce(rbind, lapply(m, function(x) x[["natural"]])),
              shifted = Reduce(rbind, lapply(m, function(x) x[["shifted"]])),
              sl_weights = lapply(m, function(x) x[["sl_weights"]]))
  return(out)
}
