setup_cv <- function(data, id, V = 10) {
  out <- origami::make_folds(data, cluster_ids = id, V = V)
  if (V > 1) {
    return(out)
  }
  out[[1]]$training_set <- out[[1]]$validation_set
  out
}

get_folded_data <- function(data, folds, index) {
  out <- list()
  out[["train"]] <- data[folds[[index]]$training_set, , drop = FALSE]
  out[["valid"]] <- data[folds[[index]]$validation_set, , drop = FALSE]
  out
}
