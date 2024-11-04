tmle <- function(x, ...) {
  UseMethod("tmle")
}

tmle.LmtpWideTask <- function(task, outcome_reg, density_ratios, bootstrap, control) {
  outcome_regs <- estimate_tmle2.LmtpWideTask(task, outcome_reg, density_ratios)

  if (!bootstrap) {
    return(outcome_regs)
  }
}

estimate_tmle2.LmtpWideTask <- function(task, outcome_reg, density_ratios) {
  tml2 <- Tmle2Wide$new(task, outcome_reg, density_ratios)
  while (tml2$t > 0) {
    tml2$fluctuate(task)
  }
  list(natural = tml2$m_eps$natural,
       shifted = tml2$m_eps$shifted)
}

Tmle2Wide <- R6Class("Tmle2Wide",
  public = list(
    outcome_reg = NULL,
    density_ratios = NULL,
    m_eps = NULL,
    t = NULL,
    initialize = function(task, outcome_reg, density_ratios) {
      self$t <- task$tau

      self$outcome_reg <- outcome_reg
      self$density_ratios <- accumulate(density_ratios)

      self$m_eps <- lapply(1:2, \(x) matrix(0, nrow = task$nrow(), ncol = task$tau + 1))
      names(self$m_eps) <- c("natural", "shifted")
      self$m_eps$shifted[, self$t + 1] <- task$data()[[last(task$col_roles$Y)]]
    },
    fluctuate = function(task) {
      on.exit(task$reset())

      target <- self$m_eps$shifted[, self$t + 1]

      # Subset active rows/cols to observed at t and at risk observations
      task$obs(self$t)$at_risk(task$t)

      # Make weights
      case_weights <- task$weights()
      weights <- self$density_ratios[task$active_rows, self$t] * case_weights

      intercept <- self$outcome_reg$natural[task$active_rows, self$t]

      fit <- sw(glm(target ~ offset(qlogis(intercept)), weights = weights, family = "binomial"))

      # Subset active rows/cols to observed at t-1 and at risk observations
      task$reset()
      task$obs(self$t - 1)$at_risk(self$t)

      # Update natural
      self$m_eps$natural[task$active_rows, self$t] <-
        bound(plogis(qlogis(self$outcome_reg$natural[task$active_rows, self$t]) + coef(fit)))

      # Update shifted
      self$m_eps$shifted[task$active_rows, self$t] <-
        bound(plogis(qlogis(self$outcome_reg$shifted[task$active_rows, self$t]) + coef(fit)))

      # iterate back through time
      self$iterate()
    },
    iterate = function() {
      self$t <- self$t - 1
      invisible(self)
    }
  )
)

crossfit_tmle <- function(x, ...) {
  UseMethod("crossfit_tmle")
}

#' @export
crossfit_tmle.LmtpWideTask <- function(task, density_ratios, learners, control, pb) {
  ans <- vector("list", length = task$nfolds())

  # Calculate cumulative density ratios
  ratios <- accumulate(density_ratios)

  for (fold in seq_along(task$folds)) {
    train <- task$training(fold)
    valid <- task$validation(fold)
    ratios_train <- ratios[task$folds[[fold]]$training_set, , drop = FALSE]

    ans[[fold]] <- future::future(
      estimate_tmle.LmtpWideTask(train, valid, ratios_train, learners, control, pb),
      seed = TRUE
    )
  }

  ans <- future::value(ans)

  list(natural = recombine(rbind_depth(ans, "natural"), task$folds),
       shifted = recombine(rbind_depth(ans, "shifted"), task$folds),
       fits = lapply(ans, \(x) x[["fits"]]))
}

estimate_tmle.LmtpWideTask <- function(train, valid, density_ratios, learners, control, pb) {
  tml <- TmleWide$new(train, valid, density_ratios)
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

    tml$add_m(fit, train)
    tml$add_m(fit, valid)
    pseudo <- tml$fluctuate(train, valid)

    train$modify(target, pseudo[train$active_rows])

    pb()
  }

  list(natural = tml$m_valid$natural,
       shifted = tml$m_valid$shifted,
       fits = fits)
}

TmleWide <- R6Class("TmleWide",
  inherit = EstimatorWide,
  public = list(
    initialize = function(train, valid, density_ratios) {
      super$initialize(train, valid, density_ratios)
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
      self$iterate()

      pseudo
    }
  )
)
