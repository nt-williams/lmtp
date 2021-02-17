setup_cv <- function(data, id, V = 10) {
  out <- origami::make_folds(data, cluster_ids = id, V = V)
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

get_folded_weights <- function(weights, folds) {
  out <- list()
  for (i in 1:length(folds)) {
    out[[i]] <- weights[folds[[i]]$training_set]
  }
  return(out)
}

cf_r <- function(data, shift, V, trt, cens, deterministic, tau,
                 node_list, learners, pb, weights_r, SL_folds) {
  fopts <- options("lmtp.bound", "lmtp.trt.length")
  out <- list()
  for (i in 1:V) {
    out[[i]] <- future::future({
      options(fopts)
      estimate_r(data[[i]]$train, data[[i]]$valid, trt, cens, deterministic,
                 shift, tau, node_list, learners, pb, weights_r[[i]], SL_folds)
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }
  out <- future::value(out)
  return(out)
}

cf_sub <- function(data, shifted, V, outcome, node_list, C, deterministic, tau,
                   outcome_type, learners, m, pb, weights_m, SL_folds) {
  fopts <- options("lmtp.bound")
  out <- list()
  for (i in 1:V) {
    out[[i]] <- future::future({
      options(fopts)
      estimate_sub(data[[i]]$train, shifted[[i]]$train, data[[i]]$valid, shifted[[i]]$valid,
                   outcome, node_list, C, deterministic, tau, outcome_type,
                   learners, m[[i]]$valid, pb, weights_m[[i]], SL_folds)
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }
  out <- future::value(out)
  out <- list(m = Reduce(rbind, lapply(out, function(x) x[["m"]])),
              sl_weights = lapply(out, function(x) x[["sl_weights"]]))
  return(out)
}

cf_tmle <- function(data, shifted, V, outcome, node_list, C, deterministic, tau,
                    outcome_type, m_natural, m_shifted, r, learners, pb,
                    weights, weights_m, SL_folds) {
  fopts <- options("lmtp.bound")
  m <- list()
  for (i in 1:V) {
    m[[i]] <- future::future({
      options(fopts)
      estimate_tmle(data[[i]]$train, shifted[[i]]$train, data[[i]]$valid,
                    shifted[[i]]$valid, outcome, node_list, C, deterministic, tau,
                    outcome_type, m_natural[[i]], m_shifted[[i]], r[[i]],
                    learners, pb, weights[[i]], weights_m[[i]], SL_folds)
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }
  m <- future::value(m)
  out <- list(natural = Reduce(rbind, lapply(m, function(x) x[["natural"]])),
              shifted = Reduce(rbind, lapply(m, function(x) x[["shifted"]])),
              sl_weights = lapply(m, function(x) x[["sl_weights"]]))
  return(out)
}

cf_sdr <- function(data, shifted, V, outcome, node_list, C, deterministic,
                   tau, outcome_type, m_natural, m_shifted, r, learners,
                   pb, weights_m, trim, SL_folds) {
  fopts <- options("lmtp.bound")
  m <- list()
  for (i in 1:V) {
    m[[i]] <- future::future({
      options(fopts)
      estimate_sdr(data[[i]]$train, shifted[[i]]$train,data[[i]]$valid,
                   shifted[[i]]$valid, outcome, node_list, C, deterministic, tau, tau,
                   outcome_type, learners, m_natural[[i]], m_shifted[[i]],
                   r[[i]], pb, weights_m[[i]], trim, SL_folds)
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }
  m <- future::value(m)
  out <- list(natural = Reduce(rbind, lapply(m, function(x) x[["natural"]])),
              shifted = Reduce(rbind, lapply(m, function(x) x[["shifted"]])),
              sl_weights = lapply(m, function(x) x[["sl_weights"]]))
  return(out)
}
