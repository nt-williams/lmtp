
# initiates an sl3 task and corresponding learner stack
initiate_ensemble <- function(outcome_type, learners = NULL) {

  # specifying learners
  sl3_profile <- learner_defaults(outcome_type = outcome_type, learners = learners)

  # building learner stack
  out <- sl3::make_learner(sl3::Lrnr_sl,
                           learners = sl3_profile$learners,
                           metalearner = sl3_profile$metalearner,
                           keep_extra = FALSE)

  # returns
  return(out)
}

# general initiator of an sl3 task
initiate_sl3_task <- function(data, Y, X, outcome_type, id = NULL, drop = FALSE) {

  # general task
  task <-
    sl3::sl3_Task$new(
      data = data,
      covariates = X,
      outcome = Y,
      outcome_type = outcome_type,
      id = id,
      drop_missing_outcome = drop
    )

  # returns
  return(task)
}

# defines default meta learner and candidates if none specified
learner_defaults <- function(outcome_type, learners) {

  # setting meta learner
  if (outcome_type %in% c("binomial", "quasibinomial", "continuous")) {
    meta <- sl3::make_learner("Lrnr_nnls",
                              convex = TRUE)
  }
  # setting candidate learners if not specified
  if (is.null(learners) & outcome_type %in% c("binomial", "quasibinomial", "continuous")) {
    learners <- sl3::make_learner_stack(sl3::Lrnr_glm, sl3::Lrnr_mean)
  } else if (!is.null(learners)) {
    learners <- learners
  }

  # returns
  out <- list(learners = learners,
              metalearner = meta)

  return(out)
}
