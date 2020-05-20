
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
    outcome_type = NULL,
    bounds = NULL,
    folds = NULL,
    weights_m = NULL,
    weights_r = NULL,
    initialize = function(data, trt, outcome, nodes, baseline, cens, k,
                          shift, outcome_type = NULL, V = 10, bounds = NULL,
                          bound = NULL) {

      # initial checks
      check_for_variables(data, trt, outcome, baseline, nodes, cens)
      check_censoring(data, cens, final_outcome(outcome))
      check_missing_data(data, trt, nodes, baseline, cens, length(nodes))
      check_scaled_conflict(data)
      check_folds(V)

      # general setup
      self$n            <- nrow(data)
      self$tau          <- length(nodes)
      self$trt          <- check_trt_length(trt, length(nodes))
      self$determ       <- check_deterministic(outcome, length(nodes))
      self$node_list    <- create_node_list(trt, nodes, baseline, k)
      self$outcome_type <- outcome_type
      self$bounds       <- y_bounds(data[[final_outcome(outcome)]], outcome_type, bounds)
      set_lmtp_options("bound", bound)

      # cross validation setup
      self$folds <- folds <- setup_cv(data, V = V)
      self$m <-
        get_folded_data(cbind(matrix(
          nrow = nrow(data), ncol = length(nodes)
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
            cens, length(nodes)
          ), folds
        )

      self$shifted_data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(shift_data(data, trt, cens, shift),
                         scale_y_continuous(data[[final_outcome(outcome)]],
                                            y_bounds(data[[final_outcome(outcome)]],
                                                     outcome_type,
                                                     bounds))),
            cens, length(nodes)
          ), folds
        )

      self$weights_m <- hold_lrnr_weights(V)
      self$weights_r <- hold_lrnr_weights(V)
    }
  )
)

prepare_r_engine <- function(data, shifted, n) {
  out    <- rbind(data, shifted)
  out$id <- rep(1:n, 2)
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
