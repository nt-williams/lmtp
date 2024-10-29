crossfit_sdr <- function(x, ...) {
  UseMethod("crossfit_sdr")
}

#' @export
crossfit_sdr.LmtpWideTask <- function(task, density_ratios, learners, control, pb) {
  ans <- vector("list", length = task$nfolds())

  for (fold in seq_along(task$folds)) {
    train <- task$training(fold)
    valid <- task$validation(fold)
    ratios_train <- density_ratios[task$folds[[fold]]$training_set, , drop = FALSE]

    ans[[fold]] <- future::future(
      estimate_sdr.LmtpWideTask(train, valid, ratios_train, learners, control, pb),
      seed = TRUE
    )
  }

  ans <- future::value(ans)

  list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
       shifted = recombine(rbind_depth(ans, "shifted"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_sdr.LmtpWideTask <- function(train, valid, density_ratios, learners, control, pb) {
  sdr <- SdrWide$new(train, valid, density_ratios)
  fits <- vector("list", length = valid$tau)
  target <- last(train$col_roles$Y)

  for (t in valid$tau:1) {
    train$reset()
    valid$reset()

    features <- train$features("Y", t)

    # Subset active rows/cols to observed at t and at risk observations
    train$obs(t)$at_risk(t)$select(c(features, target))

    if (t != valid$tau) {
      pseudo <- sdr$transform()
      train$modify(target, pseudo[train$active_rows])
    }

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

    sdr$add_m(fit, train)
    sdr$add_m(fit, valid)

    pb()
  }

  list(natural = sdr$m_valid$natural,
       shifted = sdr$m_valid$shifted,
       fits = fits)
}

SdrWide <- R6Class("SdrWide",
  inherit = EstimatorWide,
  public = list(
    tau = NULL,
    initialize = function(train, valid, density_ratios) {
      super$initialize(train, valid, density_ratios)
      self$tau <- train$tau
    },
    transform = function() {
      # iterate back through time
      self$t <- self$t - 1

      natural <- self$m_train$natural
      shifted <- self$m_train$shifted

      natural[is.na(natural)] <- -999
      shifted[is.na(shifted)] <- -999

      r <- compute_weights(self$density_ratios, self$t + 1, self$tau)
      m <- shifted[, (self$t + 2):(self$tau + 1), drop = FALSE] - natural[, (self$t + 1):self$tau, drop = FALSE]

      rowSums(r * m, na.rm = TRUE) + shifted[, self$t + 1]
    }
  )
)
