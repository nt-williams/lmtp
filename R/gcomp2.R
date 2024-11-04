crossfit_gcomp <- function(x, ...) {
  UseMethod("crossfit_gcomp")
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
