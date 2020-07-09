
Meta <- R6::R6Class(
  "Meta",
  public = list(
    data = NULL,
    shifted_data = NULL,
    trt = NULL,
    determ = NULL,
    m = NULL,
    node_list = NULL,
    n = NULL,
    tau = NULL,
    id = NULL,
    outcome_type = NULL,
    bounds = NULL,
    folds = NULL,
    weights_m = NULL,
    weights_r = NULL,
    initialize = function(data, trt, outcome, time_vary, baseline, cens, k,
                          shift, learners_trt, learners_outcome, id,
                          outcome_type = NULL, V = 10, bounds = NULL,
                          bound = NULL) {

      tau <- determine_tau(outcome, trt, cens)

      # initial checks
      check_for_variables(data, trt, outcome, baseline, time_vary, cens)
      check_censoring(data, cens, final_outcome(outcome))
      check_missing_data(data, trt, outcome, time_vary, baseline, cens, tau)
      check_scaled_conflict(data)
      check_folds(V)
      check_time_vary(time_vary)
      check_estimation_engine(learners_trt, learners_outcome)

      # general setup
      self$n            <- nrow(data)
      self$tau          <- tau
      self$trt          <- check_trt_length(trt, time_vary, cens, tau)
      self$determ       <- check_deterministic(outcome, tau)
      self$node_list    <- create_node_list(trt, tau, time_vary, baseline, k)
      self$outcome_type <- outcome_type
      self$bounds       <- y_bounds(data[[final_outcome(outcome)]], outcome_type, bounds)
      self$id           <- data$lmtp_id <- create_ids(data, id)
      set_lmtp_options("bound", bound) # global bounding option

      # cross validation setup
      self$folds <- folds <- setup_cv(data, data[["lmtp_id"]], V)
      self$m <-
        get_folded_data(cbind(matrix(
          nrow = nrow(data), ncol = tau
        ), scale_y_continuous(data[[final_outcome(outcome)]],
                             y_bounds(data[[final_outcome(outcome)]],
                                      outcome_type,
                                      bounds))),
        folds)
      self$data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(data,
                         scale_y_continuous(data[[final_outcome(outcome)]],
                                            y_bounds(data[[final_outcome(outcome)]],
                                                     outcome_type,
                                                     bounds))),
            cens, tau), folds
        )

      self$shifted_data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(shift_data(data, trt, cens, shift),
                         scale_y_continuous(data[[final_outcome(outcome)]],
                                            y_bounds(data[[final_outcome(outcome)]],
                                                     outcome_type,
                                                     bounds))),
            cens, tau), folds
        )

      self$weights_m <- hold_lrnr_weights(V)
      self$weights_r <- hold_lrnr_weights(V)
    }
  )
)

prepare_r_engine <- function(data, shifted, n) {
  out    <- rbind(data, shifted)
  out$si <- c(rep(0, n), rep(1, n))
  return(out)
}

create_r_stacks <- function(training, validation, trt, cens, shift, t, nt, nv) {
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
