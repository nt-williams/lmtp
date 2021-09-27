setup_cv <- function(data, id, V = 10) {
  origami::make_folds(data, cluster_ids = id, V = V)
}

get_folded_data <- function(data, folds) {
  out <- list()
  for (i in 1:length(folds)) {
    out[[i]] <- list()
    out[[i]][["train"]] <- data[folds[[i]]$training_set, , drop = FALSE]
    out[[i]][["valid"]] <- data[folds[[i]]$validation_set, , drop = FALSE]
  }
  out
}

get_folded_weights <- function(weights, folds) {
  out <- list()
  for (i in 1:length(folds)) {
    out[[i]] <- weights[folds[[i]]$training_set]
  }
  out
}

cf_r <- function(data, shifted, folds, trt, cens, deterministic, tau,
                 node_list, learners, pb, weights_r,
                 intervention_type, SL_folds, trim) {
  fopts <- options("lmtp.bound", "lmtp.trt.length")
  out <- list()

  for (i in seq_along(folds)) {
    out[[i]] <- future::future({
      options(fopts)

      estimate_r(
        get_folded_data(data, folds)[[i]],
        get_folded_data(shifted, folds)[[i]],
        trt, cens, deterministic,
        tau, node_list, learners,
        pb, weights_r[[i]],
        intervention_type, SL_folds
      )
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }

  trim_ratios(recombine_ratios(future::value(out), folds), trim)
}

cf_tmle <- function(data, shifted, folds, outcome, cens,
                    deterministic, tau, node_list,
                    outcome_type, m_natural, m_shifted, ratios,
                    learners, weights, weights_m, SL_folds, pb) {
  fopts <- options("lmtp.bound")
  estims <- list()

  for (i in seq_along(folds)) {
    estims[[i]] <- future::future({
      options(fopts)

      estimate_tmle(
        get_folded_data(data, folds)[[i]]$train,
        get_folded_data(shifted, folds)[[i]]$train,
        get_folded_data(data, folds)[[i]]$valid,
        get_folded_data(shifted, folds)[[i]]$valid,
        outcome,
        node_list,
        cens,
        deterministic,
        tau,
        outcome_type,
        get_folded_data(m_natural, folds)[[i]],
        get_folded_data(m_shifted, folds)[[i]],
        get_folded_data(ratios, folds)[[i]]$train,
        learners,
        pb,
        weights[[i]],
        weights_m[[i]],
        SL_folds
      )
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }

  estims <- future::value(estims)

  list(
    natural = recombine_outcome_reg(estims, "natural", folds),
    shifted = recombine_outcome_reg(estims, "shifted", folds),
    sl_weights = recombine_sl_weights(estims)
  )
}

cf_sub <- function(data, shifted, folds, outcome,
                   node_list, cens, deterministic, tau,
                   outcome_type, learners, m,
                   pb, weights_m, SL_folds) {
  fopts <- options("lmtp.bound")
  estims <- list()

  for (i in seq_along(folds)) {
    estims[[i]] <- future::future({
      options(fopts)

      estimate_sub(
        get_folded_data(data, folds)[[i]]$train,
        get_folded_data(shifted, folds)[[i]]$train,
        get_folded_data(data, folds)[[i]]$valid,
        get_folded_data(shifted, folds)[[i]]$valid,
        outcome,
        node_list,
        cens,
        deterministic,
        tau,
        outcome_type,
        learners,
        get_folded_data(m, folds)[[i]]$valid,
        pb,
        weights_m[[i]],
        SL_folds
      )
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }

  estims <- future::value(estims)

  list(
    m = recombine_outcome_reg(estims, "m", folds),
    sl_weights = recombine_sl_weights(estims)
  )
}

cf_sdr <- function(data, shifted, folds, outcome, cens,
                   deterministic, tau, node_list,
                   outcome_type, m_natural, m_shifted, ratios,
                   learners, weights, weights_m, SL_folds, pb) {
  fopts <- options("lmtp.bound")
  estims <- list()

  for (i in seq_along(folds)) {
    estims[[i]] <- future::future({
      options(fopts)

      estimate_sdr(
        get_folded_data(data, folds)[[i]]$train,
        get_folded_data(shifted, folds)[[i]]$train,
        get_folded_data(data, folds)[[i]]$valid,
        get_folded_data(shifted, folds)[[i]]$valid,
        outcome,
        node_list,
        cens,
        deterministic,
        tau,
        tau,
        outcome_type,
        learners,
        get_folded_data(m_natural, folds)[[i]],
        get_folded_data(m_shifted, folds)[[i]],
        get_folded_data(ratios, folds)[[i]]$train,
        pb,
        weights_m[[i]],
        SL_folds
      )
    },
    seed = TRUE,
    globals = structure(TRUE, add = learners))
  }

  estims <- future::value(estims)

  list(
    natural = recombine_outcome_reg(estims, "natural", folds),
    shifted = recombine_outcome_reg(estims, "shifted", folds),
    sl_weights = recombine_sl_weights(estims)
  )
}
