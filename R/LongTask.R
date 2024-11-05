LmtpLongTask <- R6Class("LmtpLongTask",
  inherit = LmtpWideTask,
  public = list(
    n = NULL,
    initialize = function(data, shifted, A, Y, W, L, C, id, weights, outcome_type, mtp, folds) {
      super$initialize(data, shifted, A, Y, W, L, C, id, weights, outcome_type, mtp, folds)

      self$task_type <- "long"
      self$n <- super$nrow()

      self$backend <- private$long("natural")
      self$shifted <- private$long("shifted")

      private$.row_copy <- self$map_active_rows(self$backend$active_row, self$backend$time)
      private$.col_copy <- names(self$backend)

      self$reset()
      invisible(self)
    },

    training = function(fold) {
      x <- subset(self$backend, active_row %in% self$folds[[fold]]$training_set)
      self$active_rows <- self$map_active_rows(x$active_row, x$time)
      LmtpLongTaskSplit$new(self, "train")
    },

    validation = function(fold) {
      x <- subset(self$backend, active_row %in% self$folds[[fold]]$validation_set)
      self$active_rows <- self$map_active_rows(x$active_row, x$time)
      LmtpLongTaskSplit$new(self, "valid")
    },

    at_risk = function() {
      if (!(self$outcome_type == "survival")) {
        return(invisible(self))
      }

      x <- self$backend
      risk <- ave(x$Y_1, x$lmtp_id, FUN = function(x) c(1, head(x, -1)))
      self$active_rows <- intersect(self$active_rows, which(risk == 1 & !is.na(risk)))
      invisible(self)
    },

    obs = function() {
      self$active_rows <- intersect(self$active_rows, which(self$backend$C_1 == 1))
      invisible(self)
    },

    map_active_rows = function(active_row, time) {
      self$n * (time - 1) + active_row
    }
  ),
  private = list(
    long = function(type) {
      vars <- lapply(1:self$tau, \(t) self$col_roles$time(t))
      to_bind <- lapply(vars, \(vars) super$data(type)[, vars])
      to_bind <- lapply(to_bind, \(data) setNames(data, self$col_roles$rename(names(data))))
      to_bind <- lapply(1:self$tau, function(t) {
        to_bind[[t]]$active_row <- self$active_rows
        to_bind[[t]]$time <- t
        to_bind[[t]]
      })
      longer <- do.call(rbind, to_bind)
      row.names(longer) <- NULL
      longer
    }
  )
)

LmtpLongTaskSplit <- R6Class("LmtpLongTaskSplit",
  inherit = LmtpLongTask,
  public = list(
    type = NULL,
    initialize = function(x, type) {
      self$type <- type
      self$n <- x$n
      self$backend <- x$data(reset = FALSE)
      self$shifted <- x$data("shifted")
      self$outcome_type <- x$outcome_type
      self$col_roles <- x$col_roles
      self$tau <- x$tau
      self$mtp <- x$mtp
      private$.row_copy <- 1:nrow(self$backend)
      private$.col_copy <- names(self$backend)
      self$active_rows <- private$.row_copy
      self$active_cols <- private$.col_copy
    }
  )
)
