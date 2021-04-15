
context("Various checks")

test_that("detects xyz", {
  df <- data.frame(xyz = 1:5)
  expect_error(check_scaled_conflict(df))
})

test_that("detects incorrect folds", {
  expect_error(check_folds(1))
})

test_that("variables dont exist", {
  a <- c("A", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  expect_error(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
                        cens, k = 1, shift = function(x) x + 0.5))
})

test_that("time_vary is a list", {
  a <- c("A1", "A2")
  nodes <- c("L1", "L2")
  cens <- c("C1", "C2")
  expect_error(lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
                        cens, k = 1, shift = function(x) x + 0.5))
})

test_that("variable length mismatch", {
  a <- c("A1")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  expect_error(
    lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
             cens, k = 1, shift = function(x) x + 0.5)
  )

  a <- c("A1", "A2")
  nodes <- list(c("L1"))
  expect_error(
    lmtp_sub(sim_cens[complete.cases(sim_cens), ], a, "Y",
             nodes, baseline = NULL, k = 1, shift = function(x) x + 0.5)
  )
})

test_that("no variation is caught", {
  x <- check_variation(rep(0.5, 10), sl3::make_learner(sl3::Lrnr_glm))
  y <- sl3::make_learner(sl3::Lrnr_mean)
  expect_equal(x$name, y$name)
})

test_that("only 0 and 1", {
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$Y <- sample(c(3, 4), nrow(sim_cens), replace = TRUE)
  expect_error(
    lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
             cens, k = 1, shift = function(x) x + 0.5)
  )
})

test_that("detect survival mismatch", {
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  expect_error(
    lmtp_sub(sim_cens, a, "Y", nodes, baseline = NULL,
             cens, k = 1, shift = function(x) x + 0.5,
             outcome_type = "survival")
  )

  a <- "trt"
  y <- paste0("Y.", 1:6)
  cens <- paste0("C.", 0:5)
  baseline <- c("W1", "W2")
  expect_error(
    lmtp_ipw(sim_point_surv, a, y, baseline, cens = cens,
             shift = static_binary_on, folds = 2)
  )
})

test_that("detect issues with supplying shifted", {
  a <- c("A1", "A2")
  tv <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  y <- "Y"
  sc <- shift_data(sim_cens, a, cens, function(data, trt) data[[trt]] + 0.5)
  sc$L1 <- 1
  expect_error(
    lmtp_tmle(sim_cens, a, y, time_vary = tv, cens = cens, shifted = sc, folds = 2)
  )
  sc <- shift_data(sim_cens, a, NULL, function(data, trt) data[[trt]] + 0.5)
  expect_error(
    lmtp_tmle(sim_cens, a, y, time_vary = tv, cens = cens, shifted = sc, folds = 2)
  )
})
