
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
    initialize = function(data, trt, outcome, nodes, baseline,
                          k, shift, outcome_type = NULL, bounds = NULL) {
      self$n <- nrow(data)
      self$tau <- length(nodes)
      self$node_list <- create_node_list(trt, nodes, baseline, k)
      self$outcome_type <- outcome_type
      self$bounds <- y_bounds(data[[outcome]], outcome_type, bounds)
      self$m <- cbind(matrix(nrow = nrow(data), ncol = length(nodes)), data[[outcome]])
      self$data <- add_scaled_y(data,
                                scale_y_continuous(data[[outcome]],
                                                   y_bounds(data[[outcome]],
                                                            outcome_type,
                                                            bounds)))
      self$shifted_data <-
        add_scaled_y(shift_data(data, trt, shift),
                     scale_y_continuous(data[[outcome]],
                                        y_bounds(data[[outcome]],
                                                 outcome_type,
                                                 bounds)))
    }
  )
)
