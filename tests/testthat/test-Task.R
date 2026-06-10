set.seed(56)
n <- 1000
W <- rnorm(n, 10, 5)
A <- 23 + 5*W + rnorm(n)
Y <- 7.2*A + 3*W + rnorm(n)
C <- rbinom(n, 1, 0.75)
dat <- data.frame(W, A, C, Y = ifelse(C == 1, Y, NA))

policy <- function(data, x) {
  (data[[x]] > 80)*(data[[x]] - 15) + (data[[x]] <= 80)*data[[x]]
}

test_that("make shifted works", {
  shifted <- dat
  shifted$A <- policy(dat, "A")
  shifted$C <- 1

  expect_equal(make_shifted(dat, "A", "C", policy, NULL), shifted)
  expect_equal(make_shifted(dat, "A", "C", NULL, shifted), shifted)
})

test_that("task creation with most basic scenario", {
  shifted <- make_shifted(dat, "A", "C", policy, NULL)
  task <- LmtpTask$new(dat, shifted, "A", "Y", NULL, "W", "C", NULL, Inf, NULL, "continuous", NULL, 10, NULL)

  expect_equal(task$id, 1:nrow(dat))
  expect_equal(task$survival, FALSE)
  expect_equal(task$.__enclos_env__$private$bounds, c(min(dat$Y, na.rm = TRUE), max(dat$Y, na.rm = TRUE)))
  expect_equal(task$rescale(task$natural$Y), dat$Y)
})

test_that("task creation with survival outcome", {
  sim_point_surv$id <- rep(1:100, each = 20)
  weights <- runif(nrow(sim_point_surv))
  shifted <- make_shifted(sim_point_surv, "trt", paste0("C.", 0:5), static_binary_on, NULL)
  task <- LmtpTask$new(sim_point_surv, shifted, "trt",
                       paste0("Y.", 1:6), NULL, c("W1", "W2"),
                       paste0("C.", 0:5), NULL, Inf, "id", "survival", NULL, 10, weights)

  expect_equal(task$survival, TRUE)
  expect_equal(task$natural$..i..lmtp_id, sim_point_surv$id)
  expect_equal(task$.__enclos_env__$private$bounds, c(0, 1))
  expect_equal(task$vars$N, paste0("Y.", 1:5))
  expect_equal(task$natural$Y.6, convert_to_surv(sim_point_surv$Y.6))
  expect_equal(mean(task$weights), 1)
  expect_equal(head(task$observed(sim_point_surv, 1)), rep(TRUE, 6))
  expect_equal(head(task$observed(sim_point_surv, 5)), c(F, F, T, F, F, T))
  expect_equal(head(task$is_outcome_free(sim_point_surv, 0)), rep(TRUE, 6))
  expect_equal(head(task$is_competing_risk_free(sim_point_surv, 0)), rep(TRUE, 6))
  expect_equal(head(task$is_outcome_free(sim_point_surv, 5)), c(F, F, F, T, F, F))
  expect_equal(head(task$is_competing_risk_free(sim_point_surv, 5)), rep(TRUE, 6))
})

test_that("task folds respect cluster ids when id= is supplied", {
  set.seed(1)
  n_clusters <- 20
  cid <- rep(seq_len(n_clusters), each = 10)
  n <- length(cid)
  d <- data.frame(W = rnorm(n), A = rnorm(n), Y = rnorm(n), cluster_id = cid)
  sh <- make_shifted(d, "A", NULL, function(data, x) data[[x]] - 1, NULL)

  task <- LmtpTask$new(
    data = d, shifted = sh,
    A = "A", Y = "Y", L = NULL, W = "W", C = NULL, D = NULL,
    k = Inf, id = "cluster_id", outcome_type = "continuous",
    bounds = NULL, folds = 5, weights = NULL
  )

  for (f in task$folds) {
    train_clusters <- unique(cid[f$training_set])
    val_clusters   <- unique(cid[f$validation_set])
    expect_length(intersect(train_clusters, val_clusters), 0)
  }
})
