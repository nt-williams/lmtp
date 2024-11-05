lmtp_gcomp <- function(data, trt, outcome, baseline = NULL, time_vary = NULL,
                     cens = NULL, shift = NULL, shifted = NULL,
                     k = Inf, outcome_type = c("binomial", "continuous", "survival"),
                     id = NULL, bounds = NULL, learners = "glm",
                     folds = 10, weights = NULL, control = lmtp_control()) {

  task <- LmtpLongTask$new(
    data = data,
    shifted = make_shifted(data, trt, cens, shift, shifted),
    A = trt,
    Y = outcome,
    W = baseline,
    L = time_vary,
    C = cens,
    id = id,
    weights = weights,
    outcome_type = match.arg(outcome_type),
    mtp = FALSE,
    folds = folds
  )

  crossfit_gcomp(task, learners, control)

}
