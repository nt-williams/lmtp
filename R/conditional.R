cf_conditional <- function(task, learners, mtp, control, pb) {
  out <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_conditional(
        get_folded_data(task$natural, task$folds, fold),
        get_folded_data(task$shifted, task$folds, fold),
        get_folded_data(task$conditional, task$folds, fold),
        task$trt,
        task$cens,
        task$risk,
        task$tau,
        task$node_list$trt,
        learners,
        pb,
        mtp,
        control
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    G = recombine_conditional(out, "G", task$folds),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_conditional <- function(natural, shifted, conditional, trt, cens, risk, tau, node_list, learners, pb, mtp, control) {
  G <- matrix(nrow = nrow(natural$valid), ncol = tau + 1)
  fits <- vector("list", length = tau)

  G[, tau + 1] <- 1

  for (t in 0:(tau - 1)) {
    if(t == 0) {
      cumulative_indicator <- apply(conditional$train[, 1:(tau + 1), drop = FALSE], 1, prod)
      G[, t + 1] <- mean(cumulative_indicator)
    }
    else {
      vars <- c(node_list[[t]], cens[[t]])
      data_train <- natural$train[c("lmtp_id", vars)]
      data_train$tmp_lmtp_conditional <- apply(conditional$train[, (t + 1):(tau + 1), drop = FALSE], 1, prod)

      if(all(data_train$tmp_lmtp_conditional == 1)) {
        G[, t] <- 1
      }
      else if(all(data_train$tmp_lmtp_conditional == 0)) {
        G[, t] <- 0
      }
      else {
        data_valid <- natural$valid[c("lmtp_id", vars)]

        fit <- run_ensemble(data_train,
                            "tmp_lmtp_conditional",
                            learners,
                            "binomial",
                            "lmtp_id",
                            control$.learners_conditional_folds)

        if (control$.return_full_fits) {
          fits[[t]] <- fit
        } else {
          fits[[t]] <- extract_sl_weights(fit)
        }

        pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
        pred <- bound(SL_predict(fit, data_valid[, c("lmtp_id", vars)]), .Machine$double.eps)

        G[, t + 1] <- pred
      }
    }

    pb()
  }

  list(G = G, fits = fits)
}
