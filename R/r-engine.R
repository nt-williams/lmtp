
#' Density Ratio Engine
#'
#' @param training Training data.
#' @param validation Validation data.
#' @param trt Name of exposure variable.
#' @param cens Names of cnesoring indicators.
#' @param C Matrices of censoring ratios.
#' @param shift Shift function.
#' @param tau Max time point.
#' @param node_list Node list created by \code{create_node_list()}.
#' @param learners An \code{sl3} learner stack.
#' @param pb Progress bar.
#'
#' @keywords internal
#' @export
estimate_r <- function(training, validation, trt, cens, C,
                       shift, tau, node_list, learners = NULL, pb) {

  # global setup
  nt <- nrow(training)
  nv <- nrow(validation)
  rt <- list(natural = matrix(nrow = nt, ncol = tau),
             shifted = matrix(nrow = nt, ncol = tau))
  rv <- list(natural = matrix(nrow = nv, ncol = tau),
             shifted = matrix(nrow = nv, ncol = tau))

  if (!is.null(shift)) {

    for (t in 1:tau) {

      # progress bar
      progress_progress_bar(pb)

      # setup
      train_stck <- prepare_r_engine(training, shift_data(training, trt[[t]], shift), nt)
      valid_stck <- prepare_r_engine(validation, shift_data(validation, trt[[t]], shift), nv)

      # create sl3 tasks for training and validation sets
      fit_task  <-
        initiate_sl3_task(
          subset(train_stck, rep(create_censoring_indicators(training, cens, tau)$i, 2)),
          "si", node_list[[t]], "binomial", "id"
        )
      pred_task <- suppressWarnings(initiate_sl3_task(valid_stck, "si", node_list[[t]], "binomial", "id"))
      ensemble  <- initiate_ensemble("binomial", learners)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # ratios
      pred            <- bound(predict_sl3(fit, fit_task), .Machine$double.eps)
      rat             <- pred / (1 - truncate(pred))
      rt$natural[, t] <- rat[train_stck$si == 0] * C$train[, t]
      rt$shifted[, t] <- rat[train_stck$si == 1] * C$train[, t]

      pred            <- bound(predict_sl3(fit, pred_task), .Machine$double.eps)
      rat             <- pred / (1 - truncate(pred))
      rv$natural[, t] <- rat[valid_stck$si == 0] * C$valid[, t]
      rv$shifted[, t] <- rat[valid_stck$si == 1] * C$valid[, t]
    }
  } else {

    for (t in 1:tau) {
      # progress bar
      progress_progress_bar(pb)

      # propensity
      rt$natural[, t] <- C$train[, t]
      rt$shifted[, t] <- C$train[, t]
      rv$natural[, t] <- C$valid[, t]
      rv$shifted[, t] <- C$valid[, t]
    }
  }

  out <- list(train = rt,
              valid = rv)

  # returns
  return(out)
}

#' Censoring Mechanism Engine
#'
#' @param data Full data set.
#' @param training Training data.
#' @param validation Validation data.
#' @param C Names of censoring nodes.
#' @param outcome Name of outcome node.
#' @param node_list Node list created by \code{create_node_list()}.
#' @param learners An \code{sl3} learner stack.
#'
#' @keywords internal
#' @export
estimate_c <- function(data, training, validation, C,
                       outcome, tau, node_list, learners) {

  # global setup
  out <- check_censoring(data, training, validation, C, outcome, tau)

  if (all(is.na(out$valid))) {
    for (t in 1:tau) {
      # setup
      fit_task  <- suppressWarnings(initiate_sl3_task(training, C[[t]], node_list[[t]], "binomial", drop = TRUE))
      pred_task <- suppressWarnings(initiate_sl3_task(validation, C[[t]], node_list[[t]], "binomial", drop = TRUE))
      ensemble  <- initiate_ensemble("binomial", learners)

      # run SL
      fit <- run_ensemble(ensemble, fit_task)

      # probability of not being censored training
      out$train[, t] <- mean(data[, C[[t]]]) / bound(predict_sl3(fit, fit_task), .Machine$double.eps)
      out$valid[, t] <- mean(data[, C[[t]]]) / bound(predict_sl3(fit, pred_task), .Machine$double.eps)
    }
  }

  # returns
  return(out)
}

ratio_dr <- function(ratios, V) {
  out <- list()
  for (i in 1:V) {
      out[[i]] <- list()
      out[[i]]$train <- check_extreme_ratio(
        matrix(t(apply(ratios[[i]]$train$natural, 1, cumprod)),
               nrow = nrow(ratios[[i]]$train$natural),
               ncol = ncol(ratios[[i]]$train$natural))
      )
      out[[i]]$valid <- check_extreme_ratio(
        matrix(t(apply(ratios[[i]]$valid$natural, 1, cumprod)),
               nrow = nrow(ratios[[i]]$valid$natural),
               ncol = ncol(ratios[[i]]$valid$natural))
      )
  }
  return(out)
}

ratio_ipw <- function(ratio) {
  out <- matrix(t(apply(ratio, 1, cumprod)), nrow = nrow(ratio), ncol = ncol(ratio))
  return(check_extreme_ratio(out))
}

ratio_sdr <- function(ratio, tau, max_tau) {
  out <- t(apply(ratio$natural[, (tau + 1):max_tau, drop = FALSE], 1, cumprod))
  if (tau == max_tau - 1) out <- t(out)
  return(check_extreme_ratio(out))
}
