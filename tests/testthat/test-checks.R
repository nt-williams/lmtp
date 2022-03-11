context("Argument checks")

test_that("Reserved variables", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$tmp_lmtp_stack_indicator <- sim_cens$Y

  expect_error(
    lmtp_sub(sim_cens, a, "tmp_lmtp_stack_indicator", time_vary = L, cens = cens),
    "Assertion on 'data' failed: 'lmtp_id', 'tmp_lmtp_stack_indicator', and 'tmp_lmtp_scaled_outcome' are reserved variable names."
  )
})

test_that("Incorrect folds", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, a, "Y", time_vary = L, cens = cens, folds = 0),
    "Assertion on 'folds' failed: Element 1 is not >= 1."
  )
})

test_that("Variables dont exist", {
  A <- c("A", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'c(trt, outcome, baseline, unlist(time_vary), cens, id)' failed: Must be a subset of {'L1','A1','C1','L2','A2','C2','Y'}, but is {'A','A2','Y','L1','L2','C1','C2'}.",
    fixed = TRUE
  )
})

test_that("Time_vary is a list", {
  A <- c("A1", "A2")
  L <- c("L1", "L2")
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'time_vary' failed: Must be of type 'list' (or 'NULL'), not 'character'.",
    fixed = TRUE
  )
})

test_that("Variable length mismatch", {
  A <- c("A1")
  L <- c("L1", "L2")
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'cens' failed: Must have length 1, but has length 2."
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
