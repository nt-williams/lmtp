#' @importFrom nnls nnls
run_ensemble <- function(data, y, learners, outcome_type, id, folds) {
  family <- ifelse(outcome_type == "binomial", binomial(), gaussian())
  cv_control <- SuperLearner::SuperLearner.CV.control(V = folds)
  features <- setdiff(names(data), c(id, y))
  X <- data[, features, drop = FALSE]
  Y <- data[[y]]
  fit <- SuperLearner::SuperLearner(
    Y, X, family = family[[1]], SL.library = learners,
    id = data[[id]], method = "method.NNLS",
    env = environment(SuperLearner::SuperLearner),
    cvControl = cv_control
  )

  class(fit) <- append("lmtp_ensemble", class(fit))
  fit
}

#' @export
predict.lmtp_ensemble <- function(object, newdata, tol = .Machine$double.eps, ...) {
  pred <- NextMethod("predict", newdata = newdata[, object$varNames], onlySL = TRUE)
  pred <- pred$pred[, 1]
  if (is.null(tol)) {
    return(pred)
  }
  bound(pred, tol)
}
