recombine <- function(x, ...) {
  UseMethod("recombine")
}

#' @export
recombine.matrix <- function(x, folds, ...) {
  i <- Reduce(c, lapply(folds, function(x) x[["validation_set"]]))
  x[order(i), , drop = FALSE]
}

rbind_depth <- function(x, n) {
  Reduce(rbind, lapply(x, function(x) x[[n]]))
}

rbind_depth_2 <- function(x, n, t) {
  vals <- lapply(x, function(x) x[[n]])
  Reduce("rbind", lapply(vals, function(x) x[[t]]))
}
