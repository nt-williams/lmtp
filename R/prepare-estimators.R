
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
    initialize = function(data, trt, outcome, nodes, baseline, cens,
                          k, shift, outcome_type = NULL, bounds = NULL) {

      check_scaled_conflict(data)

      self$n <- nrow(data)
      self$tau <- length(nodes)
      self$node_list <- create_node_list(trt, nodes, baseline, k)
      self$outcome_type <- outcome_type
      self$bounds <- y_bounds(data[[outcome]], outcome_type, bounds)
      self$m <- cbind(matrix(nrow = nrow(data), ncol = length(nodes)), data[[outcome]])
      self$data <-
        fix_censoring_ind(
          add_scaled_y(data,
                       scale_y_continuous(data[[outcome]],
                                          y_bounds(data[[outcome]],
                                                   outcome_type,
                                                   bounds))),
          cens, length(nodes)
        )
      self$shifted_data <-
        fix_censoring_ind(
          add_scaled_y(shift_data(data, trt, shift),
                       scale_y_continuous(data[[outcome]],
                                          y_bounds(data[[outcome]],
                                                   outcome_type,
                                                   bounds))),
          cens, length(nodes)
        )
    }
  )
)

prepare_r_engine <- function(data, shifted, n) {
  out    <- rbind(data, shifted)
  out$id <- rep(1:n, 2)
  out$si <- c(rep(0, n), rep(1, n))
  return(out)
}
