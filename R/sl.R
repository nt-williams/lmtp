check_variation <- function(outcome, learners) {
  if (sd(outcome) < .Machine$double.eps) {
    return("SL.mean")
  }
  learners
}

run_ensemble <- function(data, y, learners, outcome_type, id, metalearner, folds) {
  fit <- mlr3superlearner(data = data,
                          target = y,
                          library = learners,
                          metalearner = metalearner,
                          outcome_type = outcome_type,
                          folds = folds,
                          group = id)

  # if (all(fit$weights$coef == 0)) {
  #   warning("SuperLearner fit failed. Trying main-effects GLM.", call. = FALSE)
  #   tmp <- data[, setdiff(names(data), c(y, id))]
  #   tmp$lmtp_tmp_outcome_vector <- data[[y]]
  #   fit <- glm(lmtp_tmp_outcome_vector ~ ., data = tmp, family = ifelse(outcome_type == "continuous", "gaussian", "binomial"))
  # }
  fit
}

SL_predict <- function(fit, newdata) {
  if (inherits(fit, "glm")) {
    return(as.vector(predict(fit, newdata, type = "response")))
  }
  predict(fit, newdata)
}
