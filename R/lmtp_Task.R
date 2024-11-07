#' @importFrom R6 R6Class
LmtpTask <- R6::R6Class(
  "LmtpTask",
  public = list(
    natural = NULL,
    shifted = NULL,
    vars = NULL,
    id = NULL,
    n = NULL,
    tau = NULL,
    outcome_type = NULL,
    survival = NULL,
    folds = NULL,
    weights = NULL,
    initialize = function(data, shifted, A, Y, L, W, C, k,
                          id, outcome_type, V, weights) {
      # Identify tau
      self$tau <- private$tau_is(Y, A)
      self$n <- nrow(data)

      # Create Vars object
      self$vars <- LmtpVars$new(W, L, A, C, Y, self$tau, k)

      # Set outcome types
      self$outcome_type <- ifelse(outcome_type %in% c("binomial", "survival"), "binomial", "continuous")
      self$survival <- outcome_type == "survival"

      # Set outcome bounds
      private$bounds <- private$y_bounds(data[[final_outcome(Y)]])

      # Add cluster IDs
      self$id <- private$add_ids(data, id)

      # Modify data.frames to work with estimators
      self$natural <- private$as_lmtp_data(data)
      self$shifted <- private$as_lmtp_data(shifted)

      # Make folds for cross-fitting
      self$folds <- private$make_folds(V)

      # Add or normalize weights
      self$weights <- private$make_weights(weights)
    },

    scale = function(x) {
      (x - private$bounds[1]) / (private$bounds[2] - private$bounds[1])
    },

    rescale = function(x) {
      (x*(private$bounds[2] - private$bounds[1])) + private$bounds[1]
    }
  ),
  private = list(
    bounds = NULL,
    tau_is = function(Y, A) {
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
        folds <- origami::make_folds(self$natural, V = V, strata_ids = strata)
      } else {
        folds <- origami::make_folds(self$natural, cluster_ids = id, V = V)
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
      c(min(x, na.rm = T), max(x, na.rm = T))
    },

    as_lmtp_data = function(x) {
      data <- data.table::copy(as.data.frame(x))
      data$._lmtp_id <- self$id
      data <- fix_censoring_ind(data, self$vars$C)

      if (self$survival) {
        for (y in self$vars$Y) {
          data.table::set(data, j = y, value = convert_to_surv(data[[y]]))
        }
      }

      Y_tau <- final_outcome(self$vars$Y)
      data[[Y_tau]] <- self$scale(data[[Y_tau]])
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
