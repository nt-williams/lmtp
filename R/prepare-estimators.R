Meta <- R6::R6Class(
  "Meta",
  public = list(
    data = NULL,
    shifted_data = NULL,
    trt = NULL,
    risk = NULL,
    m = NULL,
    node_list = NULL,
    n = NULL,
    tau = NULL,
    id = NULL,
    outcome_type = NULL,
    survival = NULL,
    bounds = NULL,
    folds = NULL,
    weights = NULL,
    weights_m = NULL,
    weights_r = NULL,
    initialize = function(data, trt, outcome, time_vary, baseline, cens, k,
                          shift, shifted, learners_trt, learners_outcome, id,
                          outcome_type = NULL, V = 10, weights = NULL,
                          bounds = NULL, bound = NULL) {
      self$tau <- determine_tau(outcome, trt, cens)

      data <- as.data.frame(fix_censoring_ind(data, cens, self$tau))

      check_for_variables(data, trt, outcome, baseline, time_vary, cens)
      check_factors(data, trt, baseline, time_vary) # for now will just be warning, custom learners can avoid this issue.
      check_censoring(data, cens, final_outcome(outcome))
      check_missing_data(data, trt, outcome, time_vary, baseline, cens, self$tau)
      check_mult_outcomes(outcome, outcome_type)
      check_is_binary(data, outcome, outcome_type)
      check_scaled_conflict(data)
      check_time_vary(time_vary)
      check_folds(V)

      self$n <- nrow(data)
      self$trt <- check_trt_length(trt, time_vary, cens, self$tau)
      self$risk <- check_at_risk(outcome, self$tau)
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
        else if (!is.null(shifted) && !is.null(data))
          check_shifted(data, shifted, outcome, baseline, time_vary, cens)
        else
          check_shifted(data, shifted, outcome, baseline, time_vary, cens)
      }

      self$m <- cbind(
        matrix(nrow = nrow(data), ncol = self$tau),
        scale_y_continuous(
          data[[final_outcome(outcome)]],
          y_bounds(
            data[[final_outcome(outcome)]],
            self$outcome_type,
            bounds
          )
        )
      )

      self$data <- fix_censoring_ind(
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

      self$shifted_data <- fix_censoring_ind(
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
        self$weights <- get_folded_weights(weights, self$folds)
      }

      self$weights_m <- hold_lrnr_weights(V)
      self$weights_r <- hold_lrnr_weights(V)
    }
  )
)

prepare_r_engine <- function(data, shifted) {
  n <- nrow(data)
  out <- rbind(data, shifted)
  out$si <- c(rep(0, n), rep(1, n))
  out
}

create_r_stacks <- function(natural, shifted, trt, cens, tau) {
  trt_tau <- trt[tau]
  use_shifted_train <- natural$train
  use_shifted_valid <- natural$valid

  if (getOption("lmtp.trt.length") == "standard" || tau == 1) {
    use_shifted_train[[trt_tau]] <- shifted$train[[trt_tau]]
    use_shifted_valid[[trt_tau]] <- shifted$valid[[trt_tau]]
  }

  if (!is.null(cens)) {
    cens_tau <- cens[tau]
    use_shifted_train[[cens_tau]] <- shifted$train[[cens_tau]]
    use_shifted_valid[[cens_tau]] <- shifted$valid[[cens_tau]]
  }

  train_stck <- prepare_r_engine(natural$train, use_shifted_train)
  valid_stck <- prepare_r_engine(natural$valid, use_shifted_valid)

  list(train = train_stck, valid = valid_stck)
}
