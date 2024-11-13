cf_curve <- function(task, ratios, learners, control, pb) {
  ans <- vector("list", length = length(task$folds))

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_curve_gcomp(task, fold, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
       shifted = recombine(rbind_depth(ans, "shifted"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_curve_gcomp <- function(task, fold, learners, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)

  # pivot into long format
  natural <- lapply(natural, \(x) pivot(x, task$vars, task$tau))
  shifted <- lapply(shifted, \(x) pivot(x, task$vars, task$tau))

  # matrix to store predictions under the intervention for the training data
  mst <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(task$natural), ncol = tau + 1))
  # matrix to store predictions under the observed for the training data
  mnt <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(task$natural), ncol = tau))
  # matrix to store predictions under the intervention for the validation data
  msv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(task$natural), ncol = tau + 1))
  # matrix to store predictions under the observed for the validation data
  mnv <- lapply(1:task$tau, \(tau) matrix(nrow = nrow(task$natural), ncol = tau))

  taus <- order(unique(natural$train$time))
  for (t in seq_along(taus)) {
    mst[[t]][, t + 1] <- subset(natural$train, time == taus[t])$..i..Y_1
    msv[[t]][, t + 1] <- subset(natural$valid, time == taus[t])$..i..Y_1
  }

  fits <- vector("list", length = task$tau)
  for (t in 1:task$tau) {
    at_risk <- natural$train$..i..N == 1
    at_risk[is.na(at_risk)] <- TRUE
    observed <- natural$train$..i..C_1 == 1
    time <- as.numeric(natural$train$time) >= t

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

    mnt <- update_m(mnt, fit, natural$train, t)
    mst <- update_m(mst, fit, shifted$train, t)
    mnv <- update_m(mnv, fit, natural$valid, t)
    msv <- update_m(mst, fit, shifted$valid, t)

    natural$train[as.numeric(shifted$train$time) >= t, "..i..Y_1"] <-
      predict_long(fit, shifted$train, t)

    pb()
  }

  list(
    natural = mnv,
    shifted = msv,
    fits = fits
  )
}

predict_long <- function(fit, newdata, t) {
  time <- as.numeric(newdata$time) >= t
  at_risk <- newdata$..i..N[time] == 1
  at_risk[is.na(at_risk)] <- TRUE
  pred <- predict(fit, newdata[time, ], 1e-05)
  pred[!at_risk] <- 0
  pred
}

update_m <- function(m, fit, newdata, t) {
  pred <- predict_long(fit, newdata, t)
  time <- as.numeric(newdata$time) >= t
  tau <- max(as.numeric(newdata$time))
  pred <- split(pred, newdata$time[time])
  for (x in t:tau) {
    m[[x]][, x - t + 1] <- pred[[as.character(x)]]
  }
  m
}
