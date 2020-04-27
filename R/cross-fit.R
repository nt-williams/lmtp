
setup_cv <- function(data, V = 10) {
  out <- origami::make_folds(data, V = V)
  return(out)
}

get_folded_data <- function(data, folds) {
  out <- list()
  for (i in 1:length(fold_index)) {
    out[[i]] <- list()
    out[[i]][["train"]] <- data[fold_index[[i]]$training_set, ]
    out[[i]][["valid"]] <- data[fold_index[[i]]$validation_set, ]
  }
  return(out)
}
