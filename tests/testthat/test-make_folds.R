context("make_folds")

n <- 20
data <- data.frame(x = seq_len(n))

# helpers
all_rows <- function(folds) sort(unlist(lapply(folds, `[[`, "validation_set")))
disjoint <- function(folds) {
  val_sets <- lapply(folds, `[[`, "validation_set")
  length(unlist(val_sets)) == length(unique(unlist(val_sets)))
}

test_that("returns a list of length V", {
  folds <- make_folds(data, V = 5)
  expect_length(folds, 5)
})

test_that("each fold has training_set and validation_set", {
  folds <- make_folds(data, V = 5)
  for (f in folds) {
    expect_named(f, c("training_set", "validation_set"))
  }
})

test_that("validation sets partition all rows", {
  folds <- make_folds(data, V = 5)
  expect_equal(all_rows(folds), seq_len(n))
})

test_that("validation sets are disjoint", {
  folds <- make_folds(data, V = 5)
  expect_true(disjoint(folds))
})

test_that("training set is complement of validation set", {
  folds <- make_folds(data, V = 5)
  for (f in folds) {
    expect_equal(sort(c(f$training_set, f$validation_set)), seq_len(n))
  }
})

test_that("V = 1 returns single fold with training = validation = all rows", {
  folds <- make_folds(data, V = 1)
  expect_length(folds, 1)
  expect_equal(folds[[1]]$training_set, seq_len(n))
  expect_equal(folds[[1]]$validation_set, seq_len(n))
})

test_that("V = n produces LOO folds (each validation set has one row)", {
  folds <- make_folds(data, V = n)
  expect_length(folds, n)
  for (f in folds) {
    expect_length(f$validation_set, 1)
    expect_length(f$training_set, n - 1)
  }
})

test_that("V > n_units warns and defaults to LOO", {
  expect_warning(
    folds <- make_folds(data, V = n + 5),
    "defaulting to leave-one-out cross-validation"
  )
  expect_length(folds, n)
  for (f in folds) {
    expect_length(f$validation_set, 1)
  }
})

test_that("strata_ids keeps each stratum spread across folds", {
  strata <- rep(c("A", "B"), each = n / 2)
  folds <- make_folds(data, V = 5, strata_ids = strata)
  for (f in folds) {
    rows <- f$validation_set
    strata_in_fold <- strata[rows]
    expect_true(any(strata_in_fold == "A"))
    expect_true(any(strata_in_fold == "B"))
  }
})

test_that("cluster_ids keeps clusters intact (no cluster split across folds)", {
  cluster_ids <- rep(seq_len(n / 2), each = 2)
  folds <- make_folds(data, V = 5, cluster_ids = cluster_ids)
  for (f in folds) {
    val_clusters <- unique(cluster_ids[f$validation_set])
    train_clusters <- unique(cluster_ids[f$training_set])
    expect_length(intersect(val_clusters, train_clusters), 0)
  }
})

test_that("cluster_ids with strata_ids: clusters intact and strata spread", {
  cluster_ids <- rep(seq_len(n / 2), each = 2)
  strata <- rep(c("A", "B"), length.out = n)

  folds <- make_folds(data, V = 4, cluster_ids = cluster_ids, strata_ids = strata)
  expect_length(folds, 4)
  expect_true(disjoint(folds))
  expect_equal(all_rows(folds), seq_len(n))
  for (f in folds) {
    val_clusters <- unique(cluster_ids[f$validation_set])
    train_clusters <- unique(cluster_ids[f$training_set])
    expect_length(intersect(val_clusters, train_clusters), 0)
  }
})

test_that("V > number of clusters warns and defaults to LOO over clusters", {
  n_clusters <- 4
  cluster_ids <- rep(seq_len(n_clusters), each = n / n_clusters)
  expect_warning(
    folds <- make_folds(data, V = n_clusters + 2, cluster_ids = cluster_ids),
    "defaulting to leave-one-out cross-validation"
  )
  expect_length(folds, n_clusters)
})
