A <- list(c("D1", "D2"))
W <- paste0("C", 1:3)
Y <- "Y"

shift <- function(data, a) {
  out <- list(
    data[[a[1]]] - 0.1,
    data[[a[2]]] - 0.5
  )
  setNames(out, a)
}

tml <- lmtp_tmle(multivariate_data, A, Y, W, shift = shift,
          outcome_type = "continuous", folds = 1,
          mtp = TRUE)

sdr <- lmtp_sdr(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1,
         mtp = TRUE)

sub <- lmtp_sub(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1)

ipw <- lmtp_ipw(multivariate_data, A, Y, W, shift = shift,
         outcome_type = "continuous", folds = 1,
         mtp = TRUE)

test_that("Multivariate intervention fidelity", {
  expect_equal(-0.5966864, tml$theta, tolerance = 0.25)
  expect_equal(-0.5966864, sdr$theta, tolerance = 0.25)
  expect_equal(-0.5966864, sub$theta, tolerance = 0.25)
  expect_equal(-0.5966864, ipw$theta, tolerance = 0.3)
})
