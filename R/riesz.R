cf_riesz <- function(task, G, control, pb) {
  out <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_riesz(get_folded_data(task$natural, task$folds, fold),
                     get_folded_data(task$shifted, task$folds, fold),
                     task$trt,
                     task$cens,
                     task$risk,
                     get_folded_data(task$conditional, task$folds, fold),
                     get_folded_data(G, task$folds, fold),
                     task$tau,
                     task$node_list$trt,
                     control,
                     pb)
    },
    seed = TRUE)
  }

  out <- future::value(out)
  recombine_ratios(out, task$folds)
}

estimate_riesz <- function(natural,
                           shifted,
                           trt,
                           cens,
                           risk,
                           conditional,
                           G,
                           tau,
                           node_list,
                           control,
                           pb) {
  weights <- matrix(0, nrow(natural$train), ncol = tau)
  riesz_valid <- matrix(data = 0, nrow = nrow(natural$valid), ncol = tau)
  fits <- vector("list", length = tau)

  for (t in 1:tau) {
    jrt <- censored(natural$train, cens, t)$j
    drt <- at_risk(natural$train, risk, t)
    irv <- censored(natural$valid, cens, t)$i
    jrv <- censored(natural$valid, cens, t)$j
    drv <- at_risk(natural$valid, risk, t)

    if (length(trt) > 1) {
      trt_t <- trt[[t]]
    } else {
      trt_t <- trt[[1]]
    }

    vars <- c(node_list[[t]], cens[[t]])

    shifted_train <- natural$train
    shifted_train[, trt_t] <- shifted$train[, trt_t]

    if (!is.null(cens)) {
      shifted_train[[cens[t]]] <- shifted$train[[cens[t]]]
    }

    if ((t - 1) == 0) {
      wts <- rep(1, nrow(natural$train))
    } else {
      wts <- weights[jrt & drt, t - 1]
    }

    ci <- conditional$train[jrt & drt, t]

    d_in <- length(vars)
    hidden <- ceiling(mean(c(d_in, 1)))
    hidden <- 20

    net <- riesznet::nn_ensemble(
      torch::nn_sequential(
        torch::nn_linear(d_in, 1),
        torch::nn_softplus()
      ),

      torch::nn_sequential(
        torch::nn_linear(d_in, hidden),
        torch::nn_relu(),
        torch::nn_dropout(0.4),
        torch::nn_linear(hidden, hidden),
        torch::nn_relu(),
        torch::nn_dropout(0.4),
        torch::nn_linear(hidden, hidden),
        torch::nn_relu(),
        torch::nn_dropout(0.4),
        torch::nn_linear(hidden, 1),
        torch::nn_softplus()
      ),
      torch::nn_sequential(
        torch::nn_linear(d_in, hidden),
        torch::nn_relu(),
        torch::nn_dropout(0.4),
        torch::nn_linear(hidden, 1),
        torch::nn_softplus()
      )
    )

    model <- riesznet::riesznet(
      data = natural$train[jrt & drt, vars, drop = FALSE],
      shifted = list(data_1 = shifted_train[jrt & drt, vars, drop = FALSE]),
      weights = wts * (ci / G$train[jrt & drt, t]),
      .f = \(data_1) data_1,
      net = net,
      epochs = control$.epochs,
      max_lr = control$.learning_rate,
      batch_size = control$.batch_size,
      weight_decay = control$.weight_decay,
      patience = control$.patience,
      verbose = TRUE
    )

    # Return the full model object or return nothing
    if (control$.return_full_fits) {
      fits[[t]] <- model
    } else {
      fits[[t]] <- model$fit$model$modules$net.meta$parameters$weight_logits |>
        torch::nnf_softmax(dim = 2) |>
        as.numeric()
    }

    weights[jrt & drt, t] <- as.numeric(predict(model, natural$train[jrt & drt, vars, drop = FALSE]))
    riesz_valid[jrv & drv, t] <- as.numeric(predict(model, natural$valid[jrv & drv, vars, drop = FALSE]))

    pb()
  }

  list(ratios = riesz_valid, fits = fits)
}
