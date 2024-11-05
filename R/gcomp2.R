crossfit_gcomp <- function(x, ...) {
  UseMethod("crossfit_gcomp")
}

#' @export
crossfit_gcomp.LmtpLongTask <- function(task, learners, control) {
  ans <- vector("list", length = task$nfolds())

  for (fold in seq_along(task$folds)) {
    train <- task$training(fold)
    valid <- task$validation(fold)

    ans[[fold]] <- future::future(
      estimate_gcomp.LmtpLongTask(train, valid, learners, control),
      seed = TRUE
    )
  }

  ans <- future::value(ans)


  # do something
}

estimate_gcomp.LmtpLongTask <- function(train, valid, learners, control) {
  gcomp <- EstimatorLong$new(train, valid, NULL)
  fits <- vector("list", length = valid$tau)
  target <- "Y_1"

  for (l in 1:valid$tau) {
    train$reset()
    valid$reset()

    dat <- train$obs()$at_risk()$data()
    features <- setdiff(names(dat), c(target, "active_row", "C_1"))

    if (l > 1) {
      active_rows <- gcomp$map_active_rows(train, dat$active_row, dat$time)
      dat[[target]] <- gcomp$m_train$natural[active_rows, gcomp$t + 1]
    }

    fit <- run_ensemble(
      dat[dat$time >= l, c(features, target)],
      target,
      learners,
      outcome_type(train, l),
      "lmtp_id",
      control$.learners_outcome_folds,
      control$.discrete,
      control$.info
    )

    gcomp$add_m(fit, train)
    gcomp$add_m(fit, valid)
    gcomp$iterate()
  }
}

#' @export
crossfit_gcomp.LmtpWideTask <- function(task, learners, control, pb) {
  ans <- vector("list", length = task$nfolds())

  for (fold in seq_along(task$folds)) {
    train <- task$training(fold)
    valid <- task$validation(fold)

    ans[[fold]] <- future::future(
      estimate_gcomp.LmtpWideTask(train, valid, learners, control, pb),
      seed = TRUE
    )
  }

  ans <- future::value(ans)

  list(
    pred = list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
                shifted = recombine(rbind_depth(ans, "shifted"), task$folds)),
    fits = lapply(ans, \(x) x[["fits"]])
  )
}

estimate_gcomp.LmtpWideTask <- function(train, valid, learners, control, pb) {
  gcomp <- EstimatorWide$new(train, valid, NULL)
  fits <- vector("list", length = valid$tau)
  target <- last(train$col_roles$Y)

  for (t in valid$tau:1) {
    train$reset()
    valid$reset()

    features <- train$features("Y", t)

    # Subset active rows/cols to observed at t and at risk observations
    train$obs(t)$at_risk(t)$select(c(features, target))

    fit <- run_ensemble(
      train$data(),
      target,
      learners,
      outcome_type(train, t),
      "lmtp_id",
      control$.learners_outcome_folds,
      control$.discrete,
      control$.info
    )

    fits[[t]] <- return_full_fits(fit, t, control)

    gcomp$add_m(fit, train)
    gcomp$add_m(fit, valid)

    pseudo <- gcomp$m_train$shifted[, gcomp$t]
    train$modify(target, pseudo)

    gcomp$iterate()
    pb()
  }

  list(natural = gcomp$m_valid$natural,
       shifted = gcomp$m_valid$shifted,
       fits = fits)
}
