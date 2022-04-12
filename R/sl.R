sl3_task <- function(data, Y, X, outcome_type, id = NULL, V = NULL) {
  sl3::sl3_Task$new(
    data = data,
    covariates = X,
    outcome = Y,
    outcome_type = outcome_type,
    id = id,
    drop_missing_outcome = drop,
    folds = if (is.null(V)) {
      V
    } else {
      origami::make_folds(cluster_ids = data[[id]], V = V)
    }
  )
}

check_variation <- function(outcome, learners) {
  if (sd(outcome) < .Machine$double.eps) {
    return(sl3::Lrnr_mean$new())
  }
  learners
}
