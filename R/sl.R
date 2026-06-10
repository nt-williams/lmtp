#' @importFrom nnls nnls
run_ensemble <- function(data, y, learners, outcome_type, id, folds) {
  n_ids <- length(unique(data[[id]]))
  if (folds > n_ids) {
    cli::cli_abort(
      "The number of Super Learner folds ({folds}) exceeds the number of unique IDs ({n_ids}) in the training data. Set {.code .learners_trt_folds} or {.code .learners_outcome_folds} in {.fn lmtp_control} to at most {n_ids}."
    )
  }
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

#' @export
summary.SuperLearner <- function(x, time = NULL, fold = NULL, level = NULL, ...) {
  values <- data.frame(risk = x$cvRisk, coef = x$coef)
  values <- data.frame(learner = rownames(values), values, row.names = NULL)
  data.table::data.table(time = time, fold = fold, level = level, values)
}
