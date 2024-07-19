riesz_superlearner_weights <- function(learners, task_valid) {
  risks <- lapply(learners, function(x) {
    x$loss(task_valid)
  })

  weights <- numeric(length(learners))
  weights[which.min(risks)] <- 1
  list(weights = weights, risk = risks)
}

#' @importFrom SuperRiesz super_riesz
run_riesz_ensemble <- function(learners, natural_train, shifted_train, conditional_train,
                               natural_valid, shifted_valid, conditional_valid, folds) {

  if(is.null(folds)) folds <- 5
  sl <- SuperRiesz::super_riesz(
    natural_train,
    list(shifted = shifted_train, weight = data.frame(weight = conditional_train)),
    library = learners,
    folds = folds,
    m = \(alpha, data) alpha(data("shifted")) * data("weight")[,1]
  )
  predictions = predict(sl, natural_valid) * mean(conditional_valid[, 1])

  list(
    predictions = predictions,
    fits = sl,
    coef = sl$weights,
    risk = sl$risk
  )
}
