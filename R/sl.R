
# initiates an sl3 task and corresponding learner stack
initiate_ensemble <- function(outcome_type, learners = NULL) {

  # specifying learners
  sl3_profile <- learner_defaults(outcome_type = outcome_type, learners = learners)

  # building learner stack
  out <- sl3::make_learner(sl3::Lrnr_sl,
                           learners = sl3_profile$learners,
                           metalearner = sl3_profile$metalearner)

  # returns
  return(out)
}

# general initiator of an sl3 task
initiate_sl3_task <- function(data, Y, X, outcome_type, id = NULL) {

  # general task
  task <-
    sl3::sl3_Task$new(
      data = data,
      covariates = X,
      outcome = Y,
      outcome_type = outcome_type,
      id = id
    )

  # returns
  return(task)
}

# defines default meta learner and candidates if none specified
learner_defaults <- function(outcome_type, learners) {

  # setting meta learner
  if (outcome_type %in% c("binomial", "quasibinomial", "continuous")) {
    meta <- sl3::make_learner("Lrnr_nnls")
  } else if (outcome_type == "density") {
    meta <- sl3::make_learner("Lrnr_solnp_density")
  }

  # setting candidate learners if not specified
  if (is.null(learners) & outcome_type %in% c("binomial", "quasibinomial", "continuous")) {
    learners <- c("Lrnr_glmnet", "Lrnr_glm", "Lrnr_mean")
    learners <- lapply(learners, sl3::make_learner)
    learners <- sl3::make_learner(sl3::Stack, learners)
  } else if (is.null(learners) & outcome_type == "density") {
    learners <- sl3::make_learner_stack(
      list("Lrnr_haldensify",
           n_bins = 5,
           grid_type = "equal_mass",
           lambda_seq = exp(seq(-1, -9, length = 300)))
    )
  } else if (!is.null(learners)) {
    learners <- learners
  }

  # returns
  out <- list(learners = learners,
              metalearner = meta)

  return(out)
}
