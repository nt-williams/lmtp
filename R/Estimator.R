EstimatorWide <- R6Class("EstimatorWide",
  public = list(
    m_train = NULL,
    m_valid = NULL,
    density_ratios = NULL,
    t = NULL,
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
    },
    iterate = function() {
      self$t <- self$t - 1
      invisible(self)
    }
  )
)
