A <- list(c("D1", "D2"))
W <- paste0("C", 1:3)
Y <- "Y"

shift <- function(data, a) {
  out <- list(
    data[[a[1]]] - 0.1,
    data[[a[2]]] - 0.5
  )
  setNames(out, a)
}

# truth is -0.603ish
lmtp_tmle(multivariate_data, A, Y, W, shift = shift,
          outcome_type = "continuous", folds = 1,
          intervention_type = "mtp",
          learners_outcome = c("SL.glm", "SL.earth"),
          learners_trt = c("SL.glm", "SL.earth"))

lmtp_sdr(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1,
         intervention_type = "mtp")

lmtp_sub(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1)

lmtp_ipw(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1,
         intervention_type = "mtp")
