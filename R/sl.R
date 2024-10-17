run_ensemble <- function(data, y, weights, learners, outcome_type, id, folds, discrete, info) {
  mlr3superlearner::mlr3superlearner(
    data = data,
    target = y,
    library = learners,
    outcome_type = outcome_type,
    wts = weights,
    folds = folds,
    group = id,
    discrete = discrete,
    info = info
  )
}
