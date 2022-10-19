context("Fidelity of estimators for time-varying treatment")

tmp <- sim_t4
a <- c("A_1", "A_2", "A_3", "A_4")
time_varying <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

for (i in a) {
  tmp[[i]] <- factor(tmp[[i]], levels = 0:5, ordered = TRUE)
}

d <- function(data, trt) {
  out <- list()
  a <- data[[trt]]
  for (i in 1:length(a)) {
    if (as.character(a[i]) %in% c("0", "1")) {
      out[[i]] <- as.character(a[i])
    } else {
      out[[i]] <- as.numeric(as.character(a[i])) - 1
    }
  }
  factor(unlist(out), levels = 0:5, ordered = TRUE)
}

truth <- 0.305

sub <- sw(lmtp_sub(tmp, a, "Y", time_vary = time_varying, shift = d, folds = 1))
ipw <- sw(lmtp_ipw(tmp, a, "Y", time_vary = time_varying, shift = d, intervention_type = "mtp", folds = 1))
tmle <- sw(lmtp_tmle(tmp, a, "Y", time_vary = time_varying, shift = d, intervention_type = "mtp", folds = 1))
sdr <- sw(lmtp_sdr(tmp, a, "Y", time_vary = time_varying, shift = d, intervention_type = "mtp", folds = 1))

test_that("time varying treatment fidelity, t = 4", {
  expect_equal(truth, sub$theta, tolerance = 0.025)
  expect_equal(truth, ipw$theta, tolerance = 0.05)
  expect_equal(truth, tmle$theta, tolerance = 0.025)
  expect_equal(truth, sdr$theta, tolerance = 0.025)
})
