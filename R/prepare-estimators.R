
prepare_mbased <- function(data, A, Y, nodes, baseline, k,
                           shift, outcome_type, bounds) {

  check_scaled_conflict(names(data))

  n      <- nrow(data)
  t      <- length(A)
  scaled <- scale_y_continuous(data[[Y]], outcome_type, bounds)

  out <- list(n = n,
              tau = t,
              m = cbind(matrix(nrow = n, ncol = t), data[, Y]),
              node_list = create_node_list(A, nodes, baseline, k),
              data = add_scaled_y(data, scaled),
              shifted_data = add_scaled_y(shift_data(data, A, shift), scaled),
              outcome_type = outcome_type,
              scale_meta = scaled)

  return(out)
}

prepare_rbased <- function(data, A, Y, nodes, baseline, k, shift) {
  out <- list(n = nrow(data),
              tau = length(A),
              node_list = create_node_list(A, nodes, baseline, k))

  return(out)
}

add_scaled_y <- function(data, scaled) {
  data$xyz <- scaled$scaled
  return(data)
}
