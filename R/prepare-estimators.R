
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
    initialize = function(data, A, Y, nodes, baseline,
                          k, shift, outcome_type = NULL, bounds = NULL) {
      self$n <- nrow(data)
      self$tau <- length(nodes)
      self$node_list <- create_node_list(A, nodes, baseline, k)
      self$outcome_type <- outcome_type
      self$bounds <- y_bounds(data[[Y]], outcome_type, bounds)
      self$m <-
        cbind(matrix(nrow = nrow(data), ncol = length(nodes)), data[[Y]])
      self$data <- add_scaled_y(data,
                                scale_y_continuous(data[[Y]],
                                                   y_bounds(data[[Y]],
                                                            outcome_type,
                                                            bounds)))
      self$shifted_data <- add_scaled_y(shift_data(data, A, shift),
                                        scale_y_continuous(data[[Y]],
                                                           y_bounds(data[[Y]],
                                                                    outcome_type,
                                                                    bounds)))
    }
  )
)
