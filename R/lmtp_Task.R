#' @importFrom R6 R6Class
lmtp_Task <- R6::R6Class(
  "lmtp_Task",
  public = list(
    natural = NULL,
    shifted = NULL,
    trt = NULL,
    cens = NULL,
    risk = NULL,
    node_list = NULL,
    n = NULL,
    tau = NULL,
    id = NULL,
    outcome_type = NULL,
    survival = NULL,
    bounds = NULL,
    folds = NULL,
    weights = NULL,
    initialize = function(data, trt, outcome, time_vary, baseline, cens, k, shift, shifted, id, outcome_type = NULL, V = 10, weights = NULL, bounds = NULL, bound = NULL) {
      self$tau <- determine_tau(outcome, trt)
      self$n <- nrow(data)
      self$trt <- trt
      self$risk <- risk_indicators(outcome)
      self$cens <- cens
      self$node_list <- create_node_list(trt, self$tau, time_vary, baseline, k)
      self$outcome_type <- ifelse(outcome_type %in% c("binomial", "survival"), "binomial", "continuous")
      self$survival <- outcome_type == "survival"
      self$bounds <- y_bounds(data[[final_outcome(outcome)]], self$outcome_type, bounds)
      data$lmtp_id <- create_ids(data, id)
      self$id <- data$lmtp_id
      self$folds <- setup_cv(data, data$lmtp_id, V)

      shifted <- {
        if (is.null(shifted) && !is.null(shift))
          shift_data(data, trt, cens, shift)
        else if (is.null(shifted) && is.null(shift))
          shift_data(data, trt, cens, shift)
        else {
          tmp <- shifted
          tmp$lmtp_id <- data$lmtp_id
          tmp
        }
      }

      data <- data.table::copy(data)
      shifted <- data.table::copy(shifted)

      data <- fix_censoring_ind(data, cens)
      shifted <- fix_censoring_ind(shifted, cens)

      if (self$survival) {
        for (outcomes in outcome) {
          data.table::set(data, j = outcomes, value = convert_to_surv(data[[outcomes]]))
          data.table::set(shifted, j = outcomes, value = convert_to_surv(shifted[[outcomes]]))
        }
      }

      data$tmp_lmtp_scaled_outcome <- scale_y(data[[final_outcome(outcome)]], self$bounds)
      shifted$tmp_lmtp_scaled_outcome <- data$tmp_lmtp_scaled_outcome

      self$natural <- data
      self$shifted <- shifted

      if (!is.null(weights)) {
        self$weights <- weights
      }
    }
  )
)
