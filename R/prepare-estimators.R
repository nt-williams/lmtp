
Meta <- R6::R6Class(
  "Meta",
  public = list(
    data = NULL,
    shifted_data = NULL,
    m = NULL,
    node_list = NULL,
    n = NULL,
    tau = NULL,
    outcome_type = NULL,
    bounds = NULL,
    folds = NULL,
    weights_m = NULL,
    weights_r = NULL,
    weights_c = NULL,
    initialize = function(data, trt, outcome, nodes, baseline, cens, k,
                          shift, outcome_type = NULL, V = 10, bounds = NULL,
                          bound = NULL, count_lrnrs_outcome, count_lrnrs_trt) {

      check_scaled_conflict(data)

      # general setup
      self$n            <- nrow(data)
      self$tau          <- length(nodes)
      self$node_list    <- create_node_list(trt, nodes, baseline, k)
      self$outcome_type <- outcome_type
      self$bounds       <- y_bounds(data[[outcome]], outcome_type, bounds)
      set_lmtp_options("bound", bound)

      # cross validation setup
      self$folds        <- folds <- setup_cv(data, V = V)
      self$m <-
        get_folded_data(cbind(matrix(
          nrow = nrow(data), ncol = length(nodes)
        ), scale_y_continuous(data[[outcome]],
                             y_bounds(data[[outcome]],
                                      outcome_type,
                                      bounds))),
        folds)
      self$data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(data,
                         scale_y_continuous(data[[outcome]],
                                            y_bounds(data[[outcome]],
                                                     outcome_type,
                                                     bounds))),
            cens, length(nodes)
          ), folds
        )

      self$shifted_data <-
        get_folded_data(
          fix_censoring_ind(
            add_scaled_y(shift_data(data, trt, shift),
                         scale_y_continuous(data[[outcome]],
                                            y_bounds(data[[outcome]],
                                                     outcome_type,
                                                     bounds))),
            cens, length(nodes)
          ), folds
        )

      self$weights_m <- create_lrnr_matrix(V, length(nodes), count_lrnrs_outcome)
      self$weights_r   <-
        self$weights_c <-
        create_lrnr_matrix(V, length(nodes), count_lrnrs_trt)
    }
  )
)

prepare_r_engine <- function(data, shifted, n) {
  out    <- rbind(data, shifted)
  out$id <- rep(1:n, 2)
  out$si <- c(rep(0, n), rep(1, n))
  return(out)
}
