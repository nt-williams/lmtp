LmtpTask <- R6Class("LmtpTask",
  public = list(
    task_type = NULL,
    backend = NULL,
    shifted = NULL,
    col_roles = NULL,
    folds = NULL,
    outcome_type = NULL,
    mtp = NULL,
    tau = NULL,

    initialize = function(task_type, outcome_type, mtp) {
      self$task_type <- task_type
      self$mtp <- mtp
      self$outcome_type <- outcome_type
    },

    history = function(var, t) {
      self$col_roles$history(var, t)
    },

    training = function(fold) {
      if (self$task_type == "wide") {
        self$active_rows <- self$folds[[fold]]$training_set
        LmtpWideTaskSplit$new(self, "train")
      }
    },

    validation = function(fold) {
      if (self$task_type == "wide") {
        self$active_rows <- self$folds[[fold]]$validation_set
        LmtpWideTaskSplit$new(self, "valid")
      }
    },

    data = function(x = c("natural", "shifted"), reset = TRUE) {
      i <- private$.active_rows
      j <- private$.active_cols

      if (reset) self$reset()

      if (match.arg(x) == "shifted") {
        return(self$shifted[i, j])
      }

      self$backend[i, j]
    },

    reset = function() {
      private$.active_rows <- private$.row_copy
      private$.active_cols <- private$.col_copy
      invisible(self)
    },

    select = function(cols) {
      assert_character(cols)
      private$.active_cols <- intersect(private$.active_cols, cols)
      invisible(self)
    },

    modify = function(col, x) {
      self$backend[self$active_rows, col] <- x
      invisible(self)
    },

    nrow = function() {
      nrow(self$backend[self$active_rows, ])
    },

    nfolds = function() {
      length(self$folds)
    }
  ),
  active = list(
    active_rows = function(rhs) {
      if (missing(rhs)) {
        return(private$.active_rows)
      }
      private$.active_rows <- rhs
    },

    active_cols = function(rhs) {
      if (missing(rhs)) {
        return(private$.active_cols)
      }
      private$.active_cols <- rhs
    }
  ),
  private = list(
    .row_copy = NULL,
    .col_copy = NULL,
    .active_rows = NULL,
    .active_cols = NULL,
    bounds = NULL,
    make_folds = function(folds) {
      assert_number(folds, lower = 1, finite = TRUE)

      if (folds == 1) {
        folded <- list(list(
          v = 1,
          training_set = 1:nrow(self$backend),
          validation_set = 1:nrow(self$backend)
        ))

        return(folded)
      }

      if (is.null(self$col_roles$id) & self$outcome_type == "binomial") {
        strata <- self$backend[[self$col_roles$Y]]
        strata[is.na(strata)] <- 2
        folded <- origami::make_folds(self$backend, V = folds, strata_ids = strata)
        return(folded)
      }

      if (!is.null(self$col_roles$id)) {
        folded <- origami::make_folds(self$backend, cluster_ids = self$backend[[self$col_roles$id]], V = folds)
        return(folded)
      }

      origami::make_folds(self$backend, V = folds)
    }
  )
)
