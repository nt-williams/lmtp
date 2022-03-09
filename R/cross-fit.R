setup_cv <- function(data, id, V = 10) {
  out <- origami::make_folds(data, cluster_ids = id, V = V)
  if (V > 1) {
    return(out)
  }
  out[[1]]$training_set <- out[[1]]$validation_set
  out
}

get_folded_data <- function(data, folds, index) {
  out <- list()
  out[["train"]] <- data[folds[[index]]$training_set, , drop = FALSE]
  out[["valid"]] <- data[folds[[index]]$validation_set, , drop = FALSE]
  out
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
    seed = TRUE)
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
    seed = TRUE)
  }

  estims <- future::value(estims)

  list(
    natural = recombine_outcome_reg(estims, "natural", folds),
    shifted = recombine_outcome_reg(estims, "shifted", folds),
    sl_weights = recombine_sl_weights(estims)
  )
}
