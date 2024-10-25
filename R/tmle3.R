crossfit_tmle <- function(x, ...) {
  UseMethod("crossfit_tmle")
}

#' @export
crossfit_tmle.LmtpWideTask <- function(task, density_ratios, learners, control, pb) {
  ans <- vector("list", length = task$nfolds())

  # Calculate cumulative density ratios
  ratios <- matrix(t(apply(density_ratios, 1, cumprod)),
                   nrow = nrow(density_ratios),
                   ncol = ncol(density_ratios))

  for (fold in seq_along(task$folds)) {
    train <- task$training(fold)
    valid <- task$validation(fold)
    dr_train <- ratios[task$folds[[fold]]$training_set, ]

    ans[[fold]] <- future::future({
      estimate_tmle.LmtpWideTask(train, valid, dr_train, learners, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
       shifted = recombine(rbind_depth(ans, "shifted"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_tmle.LmtpWideTask <- function(train, valid, density_ratios, learners, control, pb) {
  tml <- TmleWide$new(train, valid, density_ratios)
  fits <- vector("list", length = valid$tau)

  for (t in valid$tau:1) {
    train$reset()
    valid$reset()

    # NOT CURRENTLY USING THE SCALED OUTCOME...
    features <- train$history("Y", t + 1)
    target <- last(train$col_roles$Y)

    if (t == valid$tau) {
      outcome_type <- train$outcome_type
    } else {
      outcome_type <- "continuous"
    }

    # Subset active rows/cols to observed at t and at risk observations
    train$obs(t)$at_risk(t)$select(c(features, target))

    fit <- run_ensemble(
      train$data(),
      target,
      learners,
      outcome_type,
      "lmtp_id",
      control$.learners_outcome_folds,
      control$.discrete,
      control$.info
    )

    tml$add_m(fit, train)
    tml$add_m(fit, valid)
    pseudo <- tml$fluctuate(train, valid)

    train$modify(target, pseudo)
  }

  list(natural = tml$m_valid$natural,
       shifted = tml$m_valid$shifted,
       fits = fits)
}

TmleWide <- R6Class("TmleWide",
  public = list(
    m_train = NULL,
    m_valid = NULL,
    density_ratios = NULL,
    t = NULL,
    initialize = function(train, valid, density_ratios) {
      self$t <- train$tau

      self$density_ratios <- density_ratios

      self$m_train <- lapply(1:2, \(x) matrix(0, nrow = train$nrow(), ncol = train$tau))
      self$m_valid <- lapply(1:2, \(x) matrix(0, nrow = valid$nrow(), ncol = valid$tau))

      names(self$m_train) <- c("natural", "shifted")
      names(self$m_valid) <- names(self$m_train)
    },
    add_m = function(fit, task) {
      task$obs(self$t - 1)$at_risk(self$t)$select(fit$x)

      predn <- predict(fit, task$data(reset = FALSE), 1e-05)
      preds <- predict(fit, task$shift(self$t), 1e-05)

      if (task$type == "train") {
        self$m_train$natural[task$active_rows, self$t] <- predn
        self$m_train$shifted[task$active_rows, self$t] <- preds
      } else {
        self$m_valid$natural[task$active_rows, self$t] <- predn
        self$m_valid$shifted[task$active_rows, self$t] <- preds
      }

      task$reset()
      invisible(self)
    },
    fluctuate = function(train, valid) {
      train$reset()
      valid$reset()

      target <- last(train$col_roles$Y)

      # Subset active rows/cols to observed at t and at risk observations
      train$obs(self$t)$at_risk(self$t)

      case_weights <- NULL
      if (!is.null(train$col_roles$weights)) {
        case_weights <- train$data(reset = FALSE)[[train$col_roles$weights]]
      }
      case_weights <- case_weights %||% rep(1, train$nrow())

      weights <- self$density_ratios[train$active_rows, self$t] * case_weights

      y <- train$data(reset = FALSE)[[target]]
      intercept <- self$m_train$natural[train$active_rows, self$t]

      fit <- sw(glm(y ~ offset(qlogis(intercept)), weights = weights, family = "binomial"))

      # Subset active rows/cols to observed at t-1 and at risk observations
      train$reset()
      train$obs(self$t - 1)$at_risk(self$t)
      valid$obs(self$t - 1)$at_risk(self$t)

      self$m_valid$natural[valid$active_rows, self$t] <-
        bound(plogis(qlogis(self$m_valid$natural[valid$active_rows, self$t]) + coef(fit)))

      self$m_valid$shifted[valid$active_rows, self$t] <-
        bound(plogis(qlogis(self$m_valid$shifted[valid$active_rows, self$t]) + coef(fit)))

      # Return a new pseudo outcome
      pseudo <- vector("numeric", length = train$nrow())
      pseudo[train$active_rows] <- bound(plogis(qlogis(self$m_train$shifted[train$active_rows, self$t]) + coef(fit)))

      # iterate back through time
      self$t <- self$t - 1

      pseudo
    }
  )
)
