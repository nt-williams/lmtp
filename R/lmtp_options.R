lmtp_control <- function(...) {
  change <- list(...)
  control <- list(.bound = 1e-5,
                  .trim = 0.999,
                  .learners_trt_metalearner = "glm",
                  .learners_outcome_metalearner = "glm",
                  .learners_outcome_folds = 10,
                  .learners_trt_folds = 10,
                  .return_full_fits = FALSE)

  if (length(change) == 0) return(control)

  for (arg in names(change)) {
    control[[arg]] <- change[[arg]]
  }

  control
}
