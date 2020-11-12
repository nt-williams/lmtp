
run_ensemble <- function(Y, X, learners, outcome_type, id) {
  family <- ifelse(outcome_type == "binomial", binomial(), gaussian())
  SuperLearner::SuperLearner(Y, X, family = family[[1]], SL.library = learners,
                             id = id, method = "method.NNLS")
}
