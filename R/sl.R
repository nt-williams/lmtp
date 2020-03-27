
initiate_ensemble <- function(data, Y, X, outcome_type, learners = NULL) {

  # general task
  task <-
    sl3::sl3_Task$new(
      data = data,
      covariates = X,
      outcome = Y,
      outcome_type = outcome_type,
    )

  # specifying learners
  if (is.null(learners)) {
    learners <- c("Lrnr_glm", "Lrnr_mean")
  }

  learners <- lapply(learners, sl3::make_learner)
  learners <- sl3::make_learner(sl3::Stack, learners)

  if (outcome_type %in% c("binomial", "quasibinomial")) {
    meta <- sl3::make_learner("Lrnr_nnls")
  } else if (outcome_type == "density") {
    meta <- sl3::make_learner("Lrnr_solnp_density")
  }

  # building learner stack
  stack <- sl3::make_learner(sl3::Lrnr_sl,
                             learners = learners,
                             metalearner = meta)

  out <- list(task = task,
              stack = stack)

  return(out)
}

initiate_ensemble_prediction <- function(data, Y, X, outcome_type) {

  # general task
  task <-
    sl3::sl3_Task$new(
      data = data,
      covariates = X,
      outcome = Y,
      outcome_type = outcome_type,
    )

  return(task)
}
