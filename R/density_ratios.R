crossfit_density_ratio <- function(x, ...) {
  UseMethod("crossfit_density_ratio")
}

#' @export
crossfit_density_ratio.LmtpWideTask <- function(task, learners, control, pb) {
  ans <- vector("list", length = task$nfolds())

  if (length(learners) == 1 && learners == "mean") {
    warning("Using 'mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('glm', 'cv_glmnet')`.",
            call. = FALSE)
  }

  for (fold in seq_along(task$folds)) {
    train <- task$training(fold)
    valid <- task$validation(fold)

    ans[[fold]] <- future::future({
      estimate_density_ratio.LmtpWideTask(train, valid, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  ans <- list(ratios = recombine(rbind_depth(ans, "ratios"), task$folds),
              fits = lapply(ans, \(x) x[["fits"]]))

  ans$ratios <- trim(ans$ratios, control$.trim)
  ans
}

estimate_density_ratio.LmtpWideTask <- function(train, valid, learners, control, pb) {
  density_ratios <- matrix(0, nrow = valid$nrow(), ncol = valid$tau)
  fits <- vector("list", length = valid$tau)

  for (t in 1:valid$tau) {
    train$reset()
    valid$reset()

    features <- train$features("A", t)
    target <- "tmp_lmtp_stack_indicator"

    # Subset active rows/cols to observed at t-1 and at risk observations
    train$obs(t - 1)$at_risk(t)$select(c(features, target))

    fit <- run_ensemble(
      train$stack(t),
      target,
      learners,
      "binomial",
      "lmtp_id",
      control$.learners_trt_folds,
      control$.discrete,
      control$.info
    )

    fits[[t]] <- return_full_fits(fit, t, control)

    # Subset active rows/cols to observed at t and at risk observations
    valid$obs(t)$at_risk(t)$select(features)

    density_ratios[valid$active_rows, t] <- as_density_ratio(valid, fit, t)

    pb()
  }

  list(ratios = density_ratios, fits = fits)
}

as_density_ratio <- function(task, fit, t) {
  pred <- predict(fit, task$data(reset = FALSE))
  pred <- ifelse(task$followed_rule(t) & isFALSE(task$mtp), pmax(pred, 0.5), pred)
  pred / (1 - pmin(pred, 0.999))
}
