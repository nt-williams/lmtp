
# initiates an sl3 task and corresponding learner stack
initiate_ensemble <- function(outcome_type, learners = NULL) {

  if (getOption("lmtp.engine") == "sl3") {
    # specifying learners
    sl3_profile <- learner_defaults(outcome_type = outcome_type, learners = learners)

    # building learner stack
    out <- sl3::make_learner(sl3::Lrnr_sl,
                             learners = sl3_profile$learners,
                             metalearner = sl3_profile$metalearner,
                             keep_extra = FALSE)

    # returns
    return(out)
  } else {
    return(NULL) # NEED TO THINK ABOUT WHAT THIS SHOULD RETURN
  }

}

# general initiator of an sl3 task
initiate_sl3_task <- function(data, Y, X, outcome_type, id = NULL, drop = FALSE) {

  if (getOption("lmtp.engine") == "sl3") {
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
  } else {
    call_form     <- as.formula(paste(Y, " ~ ", paste(X, collapse = " + ")))
    call_data     <- substitute(data)
    call_out_type <- check_glm_outcome(outcome_type)
    return(call("glm", call_form, family = call_out_type, data = call_data))
  }

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

run_ensemble <- function(ensemble, task, envir) {
  if (!is.null(ensemble)) {
    ensemble$train(task)
  } else {
    eval(task, envir = envir)
  }
}

predict_sl3 <- function(object, task, envir) {
  if (getOption("lmtp.engine") == "sl3") {
    return(object$predict(task))
  } else {
    call_data <- eval(substitute(task$data), envir = envir)
    call_obj  <- substitute(object)
    out       <- eval(call("predict.glm", call_obj, newdata = call_data, type = "response"),
                envir = envir)
    return(as.vector(out))
  }
}
