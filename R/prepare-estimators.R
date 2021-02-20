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
                          shift, learners_trt, learners_outcome, id,
                          outcome_type = NULL, V = 10, weights = NULL,
                          bounds = NULL, bound = NULL) {
      self$tau <- determine_tau(outcome, trt, cens)

      data <- fix_censoring_ind(data, cens, self$tau)

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

      self$m <-
        get_folded_data(
          cbind(
            matrix(
              nrow = nrow(data),
              ncol = self$tau
            ),
            scale_y_continuous(
              data[[final_outcome(outcome)]],
              y_bounds(data[[final_outcome(outcome)]],
                       self$outcome_type,
                       bounds)
            )
          ),
          self$folds
        )

      self$data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(
              data,
              scale_y_continuous(
                data[[final_outcome(outcome)]],
                y_bounds(data[[final_outcome(outcome)]],
                         self$outcome_type,
                         bounds)
              )
            ),
            cens, self$tau
          ),
          self$folds
        )

      self$shifted_data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(
              shift_data(data, trt, cens, shift),
              scale_y_continuous(
                data[[final_outcome(outcome)]],
                y_bounds(data[[final_outcome(outcome)]],
                         self$outcome_type,
                         bounds)
              )
            ),
            cens, self$tau
          ),
          self$folds
        )

      if (!is.null(weights)) {
        self$weights <- get_folded_weights(weights, self$folds)
      }

      self$weights_m <- hold_lrnr_weights(V)
      self$weights_r <- hold_lrnr_weights(V)
    }
  )
)

prepare_r_engine <- function(data, shifted, n) {
  out <- rbind(data, shifted)
  out$si <- c(rep(0, n), rep(1, n))
  return(out)
}

create_r_stacks <- function(training, validation, trt, cens, shift, t, nt, nv) {
  # need to re-think this option setting here
  if (getOption("lmtp.trt.length") == "standard" || t == 1) {
    train_stck <- prepare_r_engine(training, shift_data(training, trt[[t]], cens[[t]], shift), nt)
    valid_stck <- prepare_r_engine(validation, shift_data(validation, trt[[t]], cens[[t]], shift), nv)
  } else if (getOption("lmtp.trt.length") == "point.wise" && t > 1) {
    train_stck <- prepare_r_engine(training, shift_data(training, trt[[t]], cens[[t]], NULL), nt)
    valid_stck <- prepare_r_engine(validation, shift_data(validation, trt[[t]], cens[[t]], NULL), nv)
  }
  out <- list(train = train_stck,
              valid = valid_stck)
  return(out)
}
