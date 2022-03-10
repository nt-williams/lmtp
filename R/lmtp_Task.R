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
    initialize = function(data, trt, outcome, time_vary, baseline, cens, k,
                          shift, shifted, id,
                          outcome_type = NULL, V = 10, weights = NULL,
                          bounds = NULL, bound = NULL) {
      self$tau <- determine_tau(outcome, trt)
      data <- as.data.frame(fix_censoring_ind(data, cens, self$tau))

      check_mult_outcomes(outcome, outcome_type)
      check_is_binary(data, outcome, outcome_type)

      self$n <- nrow(data)
      self$trt <- check_trt_length(trt, time_vary, cens, self$tau)
      self$risk <- check_at_risk(outcome, self$tau)
      self$cens <- cens
      self$node_list <- create_node_list(trt, self$tau, time_vary, baseline, k)
      self$outcome_type <- ifelse(outcome_type %in% c("binomial", "survival"), "binomial", "continuous")
      self$survival <- outcome_type == "survival"
      self$bounds <- y_bounds(data[[final_outcome(outcome)]], self$outcome_type, bounds)
      data$lmtp_id <- create_ids(data, id)
      self$id <- data$lmtp_id
      self$folds <- setup_cv(data, data$lmtp_id, V)

      set_lmtp_options("bound", bound)

      if (self$survival) {
        for (outcomes in outcome) {
          data.table::set(data, j = outcomes, value = convert_to_surv(data[[outcomes]]))
        }
      }

      shd <- {
        if (is.null(shifted) && !is.null(shift))
          shift_data(data, trt, cens, shift)
        else if (is.null(shifted) && is.null(shift))
          shift_data(data, trt, cens, shift)
        else if (!is.null(shifted) && !is.null(data)) {
          tmp <- check_shifted(data, shifted, outcome, baseline, time_vary, cens, self$survival)
          tmp$lmtp_id <- create_ids(tmp, id)
          tmp
        }
        else {
          tmp <- check_shifted(data, shifted, outcome, baseline, time_vary, cens, self$survival)
          tmp$lmtp_id <- create_ids(tmp, id)
          tmp
        }
      }

      self$natural <- fix_censoring_ind(
        add_scaled_y(
          data, scale_y_continuous(
            data[[final_outcome(outcome)]],
            y_bounds(
              data[[final_outcome(outcome)]],
              self$outcome_type,
              bounds
            )
          )
        ), cens, self$tau
      )

      self$shifted <- fix_censoring_ind(
        add_scaled_y(
          shd, scale_y_continuous(
            data[[final_outcome(outcome)]],
            y_bounds(
              data[[final_outcome(outcome)]],
              self$outcome_type,
              bounds
            )
          )
        ), cens, self$tau
      )

      if (!is.null(weights)) {
        self$weights <- weights
      }
    }
  )
)
