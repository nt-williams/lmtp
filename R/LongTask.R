LmtpLongTask <- R6Class("LmtpLongTask",
  inherit = "LmtpTask",
  public = list(
    initialize = function(data, trt, outcome, covar, time,
                          cens = NULL, id = NULL,
                          outcome_type, weights = NULL, folds = 1) {
      super$initialize("long", outcome_type)

      self$col_roles <- list(
        trt = trt,
        outcome = outcome,
        covar = covar,
        time = time,
        cens = cens,
        id = id,
        weights = weights
      )

      # self$backend <- as_lmtp_long_data(data)
      self$outcome_type <- outcome_type
      self$folds <- private$make_folds(folds)

      private$.row_copy <- 1:nrow(self$backend)
      private$.col_copy <- names(self$backend)
      private$.row_roles <- private$.row_copy
      private$.col_roles <- private$.col_copy
    },

    time_geq = function(t) {

    },

    obs = function() {
      if (is.null(self$col_roles$obs)) {
        return(invisible(self))
      }

      obs <- self$backend[[self$col_roles$obs]][private$.row_roles]
      private$.row_roles <- intersect(private$.row_roles, which(obs))
      invisible(self)
    },

    at_risk = function(x = c("source", "target")) {
      pop <- self$backend[[self$col_roles$pop]][private$.row_roles]
      if (match.arg(x) == "source") {
        private$.row_roles <- intersect(private$.row_roles, which(pop == 0))
      } else {
        private$.row_roles <- intersect(private$.row_roles, which(pop == 1))
      }
      invisible(self)
    }
  )
)
