cf_curve <- function(task, ratios, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_curve_sdr(task, fold, ratios, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(natural = lapply(1:task$tau, \(t) recombine(rbind_depth_2(ans, "natural", t), task$folds)),
       shifted = lapply(1:task$tau, \(t) recombine(rbind_depth_2(ans, "shifted", t), task$folds)),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_curve_sdr <- function(task, fold, ratios, learners, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  ratios <- get_folded_data(ratios, task$folds, fold)$train

  # matrix to store predictions under the intervention for the training data
  mst <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$train), ncol = tau + 1))
  # matrix to store predictions under the observed for the training data
  mnt <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$train), ncol = tau))
  # matrix to store predictions under the intervention for the validation data
  msv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau + 1))
  # matrix to store predictions under the observed for the validation data
  mnv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(natural$valid), ncol = tau))

  # pivot into long format
  natural <- lapply(natural, \(x) pivot(x, task$vars))
  shifted <- lapply(shifted, \(x) pivot(x, task$vars))

  # Keep lagged variables as natural values for predictions
  ust <- natural$train
  As <- grep("^..i..A", names(ust), value = TRUE)
  ust[, As] <- shifted$train[, As]

  usv <- natural$valid
  usv[, As] <- shifted$valid[, As]

  taus <- order(unique(natural$train$time))
  for (t in seq_along(taus)) {
    mst[[t]][, t + 1] <- subset(natural$train, time == taus[t])$..i..Y_1
    msv[[t]][, t + 1] <- subset(natural$valid, time == taus[t])$..i..Y_1
  }

  fits <- vector("list", length = task$tau)
  for (t in 1:task$tau) {
    at_risk <- natural$train$..i..N == 1
    at_risk[is.na(at_risk)] <- FALSE
    # observed <- natural$train$..i..C_1 == 1
    observed <- !is.na(natural$train$..i..Y_1)
    time <- as.numeric(natural$train$time) <= rev(1:task$tau)[t]

    vars <- setdiff(names(natural$train), c("..i..C_1", "..i..wide_id", "..i..N"))
    if (t == task$tau) {
      vars <- setdiff(vars, "time")
    }

    fit <- run_ensemble(
      natural$train[at_risk & observed & time, vars],
      "..i..Y_1",
      learners,
      ifelse(t != 1, "continuous", task$outcome_type),
      "..i..lmtp_id",
      control$.learners_outcome_folds,
      control$.discrete,
      control$.info
    )

    if (control$.return_full_fits) {
      fits[[t]] <- fit
    } else {
      fits[[t]] <- extract_sl_weights(fit)
    }

    mnt <- update_m(mnt, fit, natural$train, t, task$tau)
    mst <- update_m(mst, fit, ust, t, task$tau)
    mnv <- update_m(mnv, fit, natural$valid, t, task$tau)
    msv <- update_m(msv, fit, usv, t, task$tau)

    if ((t + 1) <= task$tau) {
      psuedo <- unlist(lapply((t + 1):task$tau, function(x) {
        l <- x - t + 1
        tau <- ncol(mnt[[x]])
        eif(ratios, mst[[x]], mnt[[x]], l, tau)
      }))
      natural$train[, "..i..Y_1"] <- c(psuedo, rep(NA_real_, nrow(natural$train) - length(psuedo)))
    }

    pb()
  }

  list(
    natural = mnv,
    shifted = msv,
    fits = fits,
    lmtp_id = natural$valid$..i..lmtp_id
  )
}

predict_long <- function(fit, newdata, t, tau) {
  time <- as.numeric(newdata$time) <= rev(1:tau)[t]
  ans <- matrix(nrow = nrow(newdata[time, ]), ncol = 1)
  observed <- newdata$..i..C_1_lag == 1
  at_risk <- newdata$..i..N[time] == 1 & !is.na(newdata$..i..N[time])
  ans[observed[time], 1] <- predict(fit, newdata[time & observed, ], NULL)
  ans[!at_risk, 1] <- 0
  ans[, 1]
}

update_m <- function(m, fit, newdata, t, tau) {
  pred <- predict_long(fit, newdata, t, tau)
  time <- as.numeric(newdata$time) >= t

  for (l in seq_along(t:tau)) {
    x <- (t:tau)[l]
    j <- as.numeric(newdata$time[time]) == x
    m[[x]][, x - t + 1] <- pred[j]
  }
  m
}
