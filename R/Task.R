#' @importFrom R6 R6Class
LmtpTask <- R6::R6Class(
  "LmtpTask",
  public = list(
    natural = NULL,
    shifted = NULL,
    vars = NULL,
    id = NULL,
    n = NULL,
    time_horizon = NULL,
    outcome_type = NULL,
    survival = NULL,
    folds = NULL,
    weights = NULL,
    initialize = function(data, shifted, A, Y, L, W, C, D, k,
                          id, outcome_type, bounds, folds, weights) {
      # Identify the time horizon
      self$time_horizon <- private$time_horizon_is(Y, A)
      self$n <- nrow(data)

      # Create Vars object
      self$vars <- LmtpVars$new(W, L, A, C, D, Y, outcome_type, self$time_horizon, k)

      # Additional checks
      assert_numeric(weights, len = nrow(data), finite = TRUE, any.missing = FALSE, null.ok = TRUE)
      assert_number(folds, lower = 1, upper = nrow(data) - 1)

      # Set outcome types
      self$outcome_type <- ifelse(outcome_type %in% c("binomial", "survival"), "binomial", "continuous")
      self$survival <- outcome_type == "survival"

      # Set outcome bounds
      if (is.null(bounds)) {
        private$bounds <- private$y_bounds(data[, self$vars$Y])
      } else {
        private$bounds <- bounds
      }

      # Add cluster IDs
      self$id <- private$add_ids(data, id)

      # Modify data.frames to work with estimators
      self$natural <- private$as_lmtp_data(data)
      self$shifted <- private$as_lmtp_data(shifted)

      assert_lmtp_data(self)

      # Make folds for cross-fitting
      self$folds <- private$make_folds(folds)

      # Add or normalize weights
      self$weights <- private$make_weights(weights)
    },

    scale = function(x) {
      (x - private$bounds[1]) / (private$bounds[2] - private$bounds[1])
    },

    rescale = function(x) {
      (x*(private$bounds[2] - private$bounds[1])) + private$bounds[1]
    },

    observed = function(data, time) {
      if (is.null(self$vars$C)) {
        return(rep(TRUE, nrow(data)))
      }

      if (time == 0) {
        return(rep(TRUE, nrow(data)))
      }

      data[[self$vars$C[time]]] == 1
    },

    is_competing_risk_free = function(data, time) {
      if (time > self$time_horizon) {
        return(rep(TRUE, nrow(data)))
      }

      # always at risk at first time point
      if (time == 0) {
        return(rep(TRUE, nrow(data)))
      }

      if (is.null(self$vars$D)) {
        return(rep(TRUE, nrow(data)))
      }

      data[[self$vars$D[time]]] == 0
    },

    is_outcome_free = function(data, time) {
      if (time > self$time_horizon) {
        return(rep(TRUE, nrow(data)))
      }

      if (is.null(self$vars$N)) {
        return(rep(TRUE, nrow(data)))
      }

      # always at risk at first time point
      if (time == 0) {
        return(rep(TRUE, nrow(data)))
      }

      data[[self$vars$N[time]]] == 1 & !is.na(data[[self$vars$N[time]]])
    },

    is_at_risk = function(data, time) {
      if (time > self$time_horizon) {
        return(rep(TRUE, nrow(data)))
      }

      # always at risk and no competing risks
      if (is.null(self$vars$N) & is.null(self$vars$D)) {
        return(rep(TRUE, nrow(data)))
      }

      self$is_competing_risk_free(data, time - 1) & self$is_outcome_free(data, time - 1)
    }

  ),
  private = list(
    bounds = NULL,
    time_horizon_is = function(Y, A) {
      if (!(length(Y) > 1)) {
        return(length(A))
      }
      length(Y)
    },

    add_ids = function(data, id) {
      if (is.null(id)) {
        return(1:nrow(data))
      }
      data[[id]]
    },

    make_folds = function(V) {
      id <- self$natural$._lmtp_id

      if (length(unique(id)) == self$n & self$outcome_type == "binomial") {
        strata <- self$natural[[final_outcome(self$vars$Y)]]
        strata[is.na(strata)] <- 2
        folds <- make_folds(self$natural, V = V, strata_ids = strata)
      } else {
        folds <- make_folds(self$natural, cluster_ids = id, V = V)
      }

      if (V > 1) {
        return(folds)
      }

      folds[[1]]$training_set <- folds[[1]]$validation_set
      folds
    },

    y_bounds = function(x) {
      if (self$outcome_type == "binomial") {
        return(c(0, 1))
      }
      c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    },

    as_lmtp_data = function(x) {
      data <- data.table::copy(as.data.frame(x))
      data$..i..lmtp_id <- self$id
      data <- fix_censoring_ind(data, self$vars$C)

      if (self$survival) {
        for (y in c(self$vars$N, self$vars$Y)) {
          data.table::set(data, j = y, value = convert_to_surv(data[[y]]))
        }
      }

      data[, self$vars$Y] <- self$scale(data[, self$vars$Y])
      data
    },

    make_weights = function(weights) {
      if (!is.null(weights)) {
        if (is_normalized(weights)) {
          return(weights)
        } else {
          # Normalize weights
          return(weights / mean(weights))
        }
      }

      rep(1, self$n)
    }
  )
)
