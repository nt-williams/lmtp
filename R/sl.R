run_ensemble <- function(data, y, learners, outcome_type, id, folds, discrete, info) {
  mlr3superlearner::mlr3superlearner(
    data = data,
    target = y,
    library = learners,
    outcome_type = outcome_type,
    folds = folds,
    group = {
      if (length(unique(data[[id]])) == nrow(data))
        NULL
      else
        id
    },
    discrete = discrete,
    info = info
  )
}

SL_predict <- function(fit, newdata) {
  if (inherits(fit, "glm")) {
    return(as.vector(predict(fit, newdata, type = "response")))
  }
  predict(fit, newdata)
}
