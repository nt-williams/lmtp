initiate_ensemble <- function(outcome_type, learners = NULL) {
  sl3::make_learner(
    sl3::Lrnr_sl,
    learners = learners,
    metalearner = sl3::make_learner("Lrnr_nnls", convex = TRUE),
    keep_extra = FALSE
  )
}

initiate_sl3_task <- function(data, Y, X, outcome_type, id = NULL, SL_folds) {
  sl3::sl3_Task$new(
    data = data,
    covariates = X,
    outcome = Y,
    outcome_type = outcome_type,
    id = id,
    drop_missing_outcome = drop,
    folds = origami::make_folds(cluster_ids = data[[id]], V = SL_folds)
  )
}

run_ensemble <- function(ensemble, task) {
  ensemble$train(task)
}

SL_predict <- function(object, task, p = getOption("lmtp.bound")) {
  bound(object$predict(task), p)
}
