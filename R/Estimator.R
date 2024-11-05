Estimator <- R6Class("Estimator",
  public = list(
    m_train = NULL,
    m_valid = NULL,
    density_ratios = NULL,
    t = NULL,
    iterate = function() {
      self$t <- self$t - 1
      invisible(self)
    }
  )
)

EstimatorLong <- R6Class("EstimatorLong",
  inherit = Estimator,
  public = list(
    initialize = function(train, valid, density_ratios) {
      self$t <- train$tau
      self$density_ratios <- density_ratios

      self$m_train <- lapply(1:2, \(x) matrix(NA_real_, nrow = train$nrow(), ncol = self$t + 1))
      self$m_valid <- lapply(1:2, \(x) matrix(NA_real_, nrow = valid$nrow(), ncol = self$t + 1))

      names(self$m_train) <- c("natural", "shifted")
      names(self$m_valid) <- names(self$m_train)

      self$m_train$shifted[, train$tau + 1] <- train$select("Y_1")$data()
      self$m_valid$shifted[, train$tau + 1] <- valid$select("Y_1")$data()
    },
    add_m = function(fit, task) {
      browser()
      task$obs()$at_risk()

      predn <- predict(fit, task$data(reset = FALSE), 1e-05)
      preds <- predict(fit, task$data("shifted", reset = FALSE), 1e-05)

      x <- task$select(c("active_row", "time"))$data()

      active_rows <- task$map_active_rows(x$active_row, x$time)

      if (task$type == "train") {
        self$m_train$natural[active_rows, self$t] <- predn
        self$m_train$shifted[active_rows, self$t] <- preds
      } else {
        self$m_valid$natural[active_rows, self$t] <- predn
        self$m_valid$shifted[active_rows, self$t] <- preds
      }

      task$reset()
      invisible(self)
    },
    map_active_rows = function(task, active_row, time) {
      nrow(task$backend) * (time - 1) + active_row
    }
  )
)

EstimatorWide <- R6Class("EstimatorWide",
  inherit = Estimator,
  public = list(
    initialize = function(train, valid, density_ratios) {
      self$t <- train$tau
      self$density_ratios <- density_ratios

      self$m_train <- lapply(1:2, \(x) matrix(0, nrow = train$nrow(), ncol = train$tau + 1))
      self$m_valid <- lapply(1:2, \(x) matrix(0, nrow = valid$nrow(), ncol = valid$tau + 1))

      names(self$m_train) <- c("natural", "shifted")
      names(self$m_valid) <- names(self$m_train)

      self$m_train$natural[, train$tau + 1] <- train$data()[[last(train$col_roles$Y)]]
      self$m_train$shifted[, train$tau + 1] <- self$m_train$natural[, train$tau + 1]

      self$m_valid$natural[, valid$tau + 1] <- valid$data()[[last(valid$col_roles$Y)]]
      self$m_valid$shifted[, valid$tau + 1] <- self$m_valid$natural[, valid$tau + 1]
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
    }
  )
)
