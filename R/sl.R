check_variation <- function(outcome, learners) {
  if (sd(outcome) < .Machine$double.eps) {
    return("SL.mean")
  }
  learners
}

#' @importFrom nnls nnls
run_ensemble <- function(Y, X, learners, outcome_type, id, folds) {
  family <- ifelse(outcome_type == "binomial", binomial(), gaussian())
  cv_control <- SuperLearner::SuperLearner.CV.control(V = folds)
  fit <- SuperLearner::SuperLearner(
    Y, X, family = family[[1]], SL.library = learners,
    id = id, method = "method.NNLS",
    env = environment(SuperLearner::SuperLearner),
    cvControl = cv_control
  )

  if (!sum(fit$coef != 0)) {
    warning("SuperLearner fit failed. Trying main-effects GLM.", call. = FALSE)
    fit <- glm(lmtp_tmp_outcome_vector ~ ., data = cbind(lmtp_tmp_outcome_vector = Y, X), family = family[[1]])
  }
  fit
}

SL_predict <- function(fit, newdata) {
  if (inherits(fit, "glm")) {
    return(as.vector(predict(fit, newdata, type = "response")))
  }
  predict(fit, newdata)$pred[, 1]
}
