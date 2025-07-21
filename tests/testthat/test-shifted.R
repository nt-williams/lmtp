context("Fidelity of estimates with shifted data supplied")

a <- c("A1", "A2")
nodes <- list(c("L1"), c("L2"))
cens <- c("C1", "C2")
truth <- 0.88

sc <- shift_data(sim_cens, a, cens, function(data, trt) data[[trt]] + 0.5)

tmle <-
  sw(lmtp_tmle(sim_cens, a, "Y", nodes, baseline = NULL,
               cens, k = 0, shifted = sc,
               outcome_type = "binomial", folds = 2, mtp = TRUE))

sdr <-
  sw(lmtp_sdr(sim_cens, a, "Y", nodes, baseline = NULL,
              cens, k = 0, shifted = sc,
              outcome_type = "binomial", folds = 2, mtp = TRUE))

test_that("estimator fidelity with shifted data supplied", {
  expect_equal(truth, tmle$estimate@x, tolerance = 0.01)
  expect_equal(truth, sdr$estimate@x, tolerance = 0.01)
})

# Testing with multivariate exposure
A <- list(c("D1", "D2"))
W <- paste0("C", 1:3)
Y <- "Y"

shifted <- multivariate_data
shifted$D1 <- multivariate_data$D1 - 0.1
shifted$D2 <- multivariate_data$D2 - 0.5

uses_shifted <- lmtp_tmle(multivariate_data, A, Y, W, shifted = shifted,
                          outcome_type = "continuous", folds = 1, mtp = TRUE)

d <- function(data, a) {
  out = list(
    data[[a[1]]] - 0.1,
    data[[a[2]]] - 0.5
  )
  setNames(out, a)
}

uses_shift <- lmtp_tmle(multivariate_data, A, Y, W, shift = d,
                        outcome_type = "continuous", folds = 1, mtp = TRUE)

test_that("shifted works with multivariate exposure", {
  expect_true(identical(uses_shifted$estimate, uses_shift$estimate))
})
