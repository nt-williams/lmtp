cf_riesz <- function(task, module, control, pb) {
  out <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_riesz(get_folded_data(task$natural, task$folds, fold),
                     get_folded_data(task$shifted, task$folds, fold),
                     task$trt,
                     task$cens,
                     task$risk,
                     task$tau,
                     task$node_list$trt,
                     module,
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
                           tau,
                           node_list,
                           module,
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

    model <- nn_riesz(
      train = list(data = natural$train[jrt & drt, vars, drop = FALSE],
                   data_1 = shifted_train[jrt & drt, vars, drop = FALSE]),
      vars = vars,
      module = module,
      .f = \(alpha, dl) alpha(dl[["data_1"]]),
      weights = wts,
      batch_size = control$.batch_size,
      learning_rate = control$.learning_rate,
      epochs = control$.epochs,
      device = control$.device
    )

    # Return the full model object or return nothing
    if (control$.return_full_fits) {
      fits[[t]] <- model
    } else {
      fits[[t]] <- NULL
    }

    weights[jrt & drt, t] <- as.numeric(
      model(
        as_torch(
          one_hot_encode(natural$train[jrt & drt, vars, drop = FALSE]),
          device = control$.device
        )
      )
    )

    riesz_valid[jrv & drv, t] <- as.numeric(
      model(
        as_torch(
          one_hot_encode(natural$valid[jrv & drv, vars, drop = FALSE]),
          device = control$.device
        )
      )
    )

    pb()
  }

  list(ratios = riesz_valid, fits = fits)
}
