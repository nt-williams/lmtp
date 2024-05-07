run_ensemble <- function(data, y, learners, outcome_type, id, folds) {
  fit <- mlr3superlearner(data = data,
                          target = y,
                          library = learners,
                          outcome_type = outcome_type,
                          folds = folds,
                          group = id)
  fit
}

SL_predict <- function(fit, newdata) {
  if (inherits(fit, "glm")) {
    return(as.vector(predict(fit, newdata, type = "response")))
  }
  predict(fit, newdata)
}
