run_ensemble <- function(Y, X, learners, outcome_type, id, folds) {
  family <- ifelse(outcome_type == "binomial", binomial(), gaussian())
  cv_control <- SuperLearner::SuperLearner.CV.control(V = folds)
  SuperLearner::SuperLearner(
    Y, X, family = family[[1]], SL.library = learners,
    id = id, method = "method.NNLS",
    env = environment(SuperLearner::SuperLearner),
    cvControl = cv_control
  )
}

SL_predict <- function(fit, newdata, p = getOption("lmtp.bound")) {
  bound(predict(fit, newdata)$pred[, 1], p)
}
