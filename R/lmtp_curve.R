lmtp_curve <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                       cens = NULL, shift = NULL, shifted = NULL,
                       k = Inf,
                       mtp = FALSE,
                       outcome_type = c("binomial", "continuous", "survival"),
                       id = NULL,
                       learners_outcome = "glm",
                       learners_trt = "glm",
                       folds = 10,
                       weights = NULL,
                       control = lmtp_control()) {

  assert_not_data_table(data)
  assert_outcome_types(data, outcome, match.arg(outcome_type))

  # Check if the treatment is continuous and warn if MTP is false
  check_trt_type(data, unlist(trt), mtp)

  task <- LmtpTask$new(
    data = data,
    shifted = make_shifted(data, trt, cens, shift, shifted),
    A = trt,
    Y = outcome,
    L = time_vary,
    W = baseline,
    C = cens,
    k = k, id = id,
    outcome_type = match.arg(outcome_type),
    V = folds, weights = weights
  )

  # Create progress bar object
  pb <- progressr::progressor(task$tau*folds*2)

  # Estimate density ratios
  ratios <- cf_r(task, learners_trt, mtp, control, pb)

  # Estimate outcome regression
  curve <- cf_curve(task, ratios$ratios, learners_outcome, control, pb)

  theta_curve(
    task = task,
    m = list(natural = curve$natural, shifted = curve$shifted),
    r = ratios$ratios,
    fits_m = curve$fits,
    fits_r = ratios$fits,
    shift = deparse(substitute((shift)))
  )
}
