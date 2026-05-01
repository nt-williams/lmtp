make_folds <- function(data, V = 10, strata_ids = NULL, cluster_ids = NULL) {
  n <- nrow(data)
  folds <- vector("list", V)

  if (V == 1) {
    folds[[1]] <- list(
      training_set = seq_len(n), 
      validation_set = seq_len(n)
    )
    return(folds)
  }

  if (is.null(cluster_ids)) {
    unit_to_rows <- as.list(seq_len(n))
    unit_strata <- if (is.null(strata_ids)) rep("1", n) else as.character(strata_ids)
  } else {
    cluster_chr <- as.character(cluster_ids)
    unique_clusters <- unique(cluster_chr)
    unit_to_rows <- lapply(unique_clusters, function(cl) which(cluster_chr == cl))

    if (is.null(strata_ids)) {
      unit_strata <- rep("1", length(unique_clusters))
    } else {
      strata_chr <- as.character(strata_ids)
      # A cluster spanning multiple strata is assigned to its modal stratum,
      # since clusters cannot be split across folds.
      unit_strata <- vapply(unique_clusters, function(cl) {
        s <- strata_chr[cluster_chr == cl]
        tab <- table(s)
        names(tab)[which.max(tab)]
      }, character(1))
    }
  }

  n_units <- length(unit_to_rows)
  if (V > n_units) {
    warning(sprintf("`k` (%d) exceeds the number of units (%d); defaulting to leave-one-out cross-validation", V, n_units))
    V <- n_units
    folds <- vector("list", V)
  }

  fold_of_unit <- integer(n_units)
  for (s in unique(unit_strata)) {
    idx <- which(unit_strata == s)
    shuffled <- sample(idx)
    fold_of_unit[shuffled] <- rep(seq_len(V), length.out = length(shuffled))
  }

  for (j in seq_len(V)) {
    val_units <- which(fold_of_unit == j)
    train_units <- which(fold_of_unit != j)
    folds[[j]] <- list(
      training_set = sort(unlist(unit_to_rows[train_units], use.names = FALSE)),
      validation_set = sort(unlist(unit_to_rows[val_units], use.names = FALSE))
    )
  }

  folds
}
