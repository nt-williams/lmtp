A = list(c("D1", "D2"))
W = paste0("C", 1:3)
Y = "Y"

create_node_list(A, 1, NULL, W)

shift = function(data, a) {
  out = list(
    data[[a[1]]] - 0.1,
    data[[a[2]]] - 0.5
  )
  setNames(out, a)
}

# truth is -0.603ish
shift_trt(multivariate_data, A, shift)

lmtp_tmle(multivariate_data, A, Y, W, shift = shift,
          outcome_type = "continuous", folds = 1,
          intervention_type = "mtp", .SL_folds = 2)

lmtp_sdr(multivariate_data, A, Y, W, shift = shift,
          outcome_type = "continuous", folds = 1,
          intervention_type = "mtp", .SL_folds = 2)

lmtp_sub(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1, .SL_folds = 2)

lmtp_ipw(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1,
         intervention_type = "mtp", .SL_folds = 2)
