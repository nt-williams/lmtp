context("Argument checks")

test_that("'data' is a 'data.frame'", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(list(), A, "tmp_lmtp_stack_indicator", time_vary = L, cens = cens),
    "Assertion on 'data' failed: Must be of type 'data.frame', not 'list'."
  )

  sim_cens <- data.table::as.data.table(sim_cens)
  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'data' failed: Must be a 'data.frame', not a 'data.table'."
  )
})

test_that("No uncensored missing data", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$A2 <- NA_real_

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'data' failed: Missing data found in treatment and/or covariate nodes for uncensored observations."
  )
})

test_that("Reserved variables", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$tmp_lmtp_stack_indicator <- sim_cens$Y

  expect_error(
    lmtp_sub(sim_cens, A, "tmp_lmtp_stack_indicator", time_vary = L, cens = cens),
    "Assertion on 'data' failed: 'lmtp_id', 'tmp_lmtp_stack_indicator', and 'tmp_lmtp_scaled_outcome' are reserved variable names."
  )
})

test_that("Incorrect folds", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 0),
    "Assertion on 'folds' failed: Element 1 is not >= 1."
  )
})

test_that("Variables dont exist", {
  A <- c("A", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'c(trt, outcome, baseline, unlist(time_vary), cens, id)' failed: Must be a subset of {'L1','A1','C1','L2','A2','C2','Y'}, but has additional elements {'A'}.",
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

test_that("No outcome variation changes learners", {
  x <- check_variation(rep(0.5, 10), "SL.glm")
  y <- "SL.mean"
  expect_equal(x, y)
})

test_that("Only 0 and 1 in 'outcome' when binary or surival", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$Y <- sample(c(3, 4), nrow(sim_cens), replace = TRUE)
  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens),
    "Assertion on 'data' failed: Only 0 and 1 allowed in outcome variables if 'outcome_type' set to binomial or survival."
  )
})

test_that("Issues with 'outcome_type' being or not being survival", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens, outcome_type = "survival"),
    "Assertion on 'outcome' failed: Must have length >= 2, but has length 1."
  )

  A <- "trt"
  Y <- paste0("Y.", 1:6)
  cens <- paste0("C.", 0:5)
  W <- c("W1", "W2")

  expect_error(
    lmtp_ipw(sim_point_surv, A, Y, W, cens = cens, shift = static_binary_on),
    "Assertion on 'outcome' failed: Must have length 1, but has length 6."
  )
})

test_that("Issues with 'shift' function and providing 'shifted' data", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  shifted <- shift_data(sim_cens, A, cens, function(data, trt) data[[trt]] + 0.5)
  shifted$L1 <- 1

  expect_error(
    lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, shift = function(data, trt, extra) data[[trt]] + 0.5),
    "Assertion on 'shift' failed: Must have exactly 2 formal arguments, but has 3."
  )

  expect_error(
    lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, shifted = shifted),
    "Assertion on 'shifted' failed: The only columns that can be different between `data` and `shifted` are those indicated in `trt` and `cens`."
  )

  shifted <- shift_data(sim_cens, A, NULL, function(data, trt) data[[trt]] + 0.5)
  expect_error(
    lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, shifted = shifted),
    "Assertion on 'shifted' failed: Censoring variables should be 1 in 'shifted'."
  )
})

test_that("Contrast assertions", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  fit <- lmtp_sub(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 1)
  expect_error(
    lmtp_contrast(fit, ref = 0.1),
    "Assertion on 'fits' failed: Contrasts not implemented for substitution/IPW estimators."
  )

  fit <- lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 1)
  expect_error(
    lmtp_contrast(fit, ref = c(0.1, 0.2)),
    "Assertion on 'ref' failed: Must either be a single numeric value or another lmtp object."
  )

  expect_error(
    lmtp_contrast(fit, "Not lmtp object", ref = 0.1),
    "Assertion on 'fits' failed: Objects must be of type 'lmtp'."
  )

  fit <- lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 1, outcome_type = "continuous")
  expect_error(
    lmtp_contrast(fit, ref = 0.1, type = "rr"),
    "Assertion on 'type' failed: 'rr' specified but one or more outcome types are not 'binomial' or 'survival'."
  )
})
