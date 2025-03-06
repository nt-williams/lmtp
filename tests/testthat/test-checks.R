context("Argument checks")

test_that("'data' is a 'data.frame'", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    lmtp_tmle(list(), A, "tmp_lmtp_stack_indicator", time_vary = L, cens = cens),
    "Assertion on 'data' failed: Must be of type 'data.frame', not 'list'."
  )

  expect_error(
    lmtp_tmle(data.table::as.data.table(sim_cens), A, "Y", time_vary = L, cens = cens),
    "Assertion on 'data' failed: Must be a 'data.frame', not a 'data.table'."
  )
})

test_that("No uncensored missing data", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$A2 <- NA_real_

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens)),
    "Assertion on 'self' failed: Missing data found in treatment and/or covariate nodes for uncensored observations."
  )
})

test_that("Incorrect folds", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 0)),
    "Assertion on 'folds' failed: Element 1 is not >= 1."
  )
})

test_that("Variables dont exist", {
  A <- c("A", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens)),
    "Assertion on 'av' failed: Must be a subset of {'L1','A1','C1','L2','A2','C2','Y'}, but has additional elements {'A'}.",
    fixed = TRUE
  )
})

test_that("Time_vary is a list", {
  A <- c("A1", "A2")
  L <- c("L1", "L2")
  cens <- c("C1", "C2")

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens)),
    "Assertion on 'L' failed: Must be of type 'list' (or 'NULL'), not 'character'.",
    fixed = TRUE
  )
})

test_that("Variable length mismatch", {
  A <- c("A1")
  L <- c("L1", "L2")
  cens <- c("C1", "C2")

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens)),
    "Assertion on 'C' failed: Must have length 1, but has length 2."
  )
})

test_that("Only 0 and 1 in 'outcome' when binary or survival", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  sim_cens$Y <- sample(c(3, 4), nrow(sim_cens), replace = TRUE)
  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens)),
    "Assertion on 'data' failed: Only 0 and 1 allowed in outcome variables if 'outcome_type' set to binomial or survival."
  )
})

test_that("Issues with 'outcome_type' being or not being survival", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, outcome_type = "survival")),
    "Assertion on 'Y' failed: Must have length >= 2, but has length 1."
  )
})

test_that("Issues with 'shift' function and providing 'shifted' data", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  shifted <- shift_data(sim_cens, A, cens, function(data, trt) data[[trt]] + 0.5)
  shifted$L1 <- 1

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, shift = function(data, trt, extra) data[[trt]] + 0.5)),
    "Assertion on 'shift' failed: Must have exactly 2 formal arguments, but has 3."
  )

  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, shifted = shifted)),
    "Assertion on 'data' failed: The only columns that can be different between `data` and `shifted` are those indicated in `trt` and `cens`."
  )

  shifted <- shift_data(sim_cens, A, NULL, function(data, trt) data[[trt]] + 0.5)
  expect_error(
    sw(lmtp_tmle(sim_cens, A, "Y", time_vary = L, cens = cens, shifted = shifted)),
    "Assertion on 'data' failed: Censoring variables should be 1 in 'shifted'."
  )
})

test_that("Contrast assertions", {
  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")

  fit <- sw(lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 1))
  expect_error(
    lmtp_contrast(fit, ref = c(0.1, 0.2)),
    "Assertion on 'ref' failed: Must either be a single numeric value or another lmtp object."
  )

  expect_error(
    lmtp_contrast(fit, "Not lmtp object", ref = 0.1),
    "Assertion on 'fits' failed: Objects must be of type 'lmtp'."
  )

  fit <- sw(lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = cens, folds = 1, outcome_type = "continuous"))
  expect_error(
    lmtp_contrast(fit, ref = 0.1, type = "rr"),
    "Assertion on 'type' failed: 'rr' specified but one or more outcome types are not 'binomial' or 'survival'."
  )
})

test_that("'bounds' issues", {
  set.seed(56)
  n <- 1000
  W <- rnorm(n, 10, 5)
  A <- 23 + 5*W + rnorm(n)
  Y <- 7.2*A + 3*W + rnorm(n)
  ex1_dat <- data.frame(W, A, Y)

  expect_error(
    lmtp_tmle(ex1_dat, "A", "Y", "W", shift = NULL,
              outcome_type = "continuous", folds = 1, mtp = TRUE,
              bounds = c(1100, -50)),
    "Assertion on 'bounds' failed: Must be sorted."
  )

  expect_error(
    lmtp_tmle(ex1_dat, "A", "Y", "W", shift = NULL,
              outcome_type = "continuous", folds = 1, mtp = TRUE,
              bounds = c(-Inf, 1100)),
    "Assertion on 'bounds' failed: Must be finite."
  )

  expect_error(
    lmtp_tmle(ex1_dat, "A", "Y", "W", shift = NULL,
              outcome_type = "continuous", folds = 1, mtp = TRUE,
              bounds = c(-50, 1100, 20000)),
    "Assertion on 'bounds' failed: Must have length 2, but has length 3."
  )
})
