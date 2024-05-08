rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
set.seed(6324)
n <- 1e4
L1 <- runif(n)
A1 <- rexpit(L1) + rexpit(L1)
Y <- rnorm(n, A1, 0.05)
sim <- data.frame(L1, A1, Y)
a <- c("A1")
time_vary <- list(c("L1"))

make_conditional <- function(a) {
  conditional0 <- matrix(
    A1 == a,
    ncol = 1
  )
}

shift <- function(data, trt) {
  ifelse(data[[trt]] == 2, 2, data[[trt]] + 1)
}

tml0 <- sw(lmtp_tmle(sim, trt = a, outcome = "Y", time_vary = time_vary,
                     conditional = make_conditional(0), outcome_type = "continuous",
                     shift = shift, folds = 5, riesz = TRUE,
                     learners_trt = list(list("nn", epochs = 100)),
                     control = lmtp_control(.learners_trt_folds = 1)))

sub0 <- sw(lmtp_sub(sim, trt = a, outcome = "Y", time_vary = time_vary,
                    conditional = make_conditional(0), outcome_type = "continuous",
                    shift = shift, folds = 1))

ipw0 <- sw(lmtp_ipw(sim, trt = a, outcome = "Y", time_vary = time_vary,
                    conditional = make_conditional(0), outcome_type = "continuous",
                    shift = shift, folds = 1, riesz = TRUE,
                    learners = list(list("nn")),
                    control = lmtp_control(.learners_trt_folds = 1)))

tml1 <- sw(lmtp_tmle(sim, trt = a, outcome = "Y", time_vary = time_vary,
                     conditional = make_conditional(1), outcome_type = "continuous",
                     shift = shift, folds = 1, riesz = TRUE,
                     learners_trt = list(list("nn", epochs = 100)),
                     control = lmtp_control(.learners_trt_folds = 1)))

sub1 <- sw(lmtp_sub(sim, trt = a, outcome = "Y", time_vary = time_vary,
                    conditional = make_conditional(1), outcome_type = "continuous",
                    shift = shift, folds = 1))

ipw1 <- sw(lmtp_ipw(sim, trt = a, outcome = "Y", time_vary = time_vary,
                    conditional = make_conditional(1), outcome_type = "continuous",
                    shift = shift, folds = 1, riesz = TRUE,
                    learners = list(list("nn", epochs = 100)),
                    control = lmtp_control(.learners_trt_folds = 1)))


tml2 <- sw(lmtp_tmle(sim, trt = a, outcome = "Y", time_vary = time_vary,
                     conditional = make_conditional(2), outcome_type = "continuous",
                     shift = shift, folds = 1, riesz = TRUE,
                     learners_trt = list(list("nn", epochs = 100)),
                     control = lmtp_control(.learners_trt_folds = 1)))

sub2 <- sw(lmtp_sub(sim, trt = a, outcome = "Y", time_vary = time_vary,
                    conditional = make_conditional(2), outcome_type = "continuous",
                    shift = shift, folds = 1))

ipw2 <- sw(lmtp_ipw(sim, trt = a, outcome = "Y", time_vary = time_vary,
                    conditional = make_conditional(2), outcome_type = "continuous",
                    shift = shift, folds = 1, riesz = TRUE,
                    learners = list(list("nn", epochs = 100)),
                    control = lmtp_control(.learners_trt_folds = 1)))


test_that("Conditional fidelity", {
  expect_equal(1, tml0$theta, tolerance = 0.1)
  expect_equal(1, sub0$theta, tolerance = 0.1)
  expect_equal(1, ipw0$theta, tolerance = 1)

  expect_equal(2, tml1$theta, tolerance = 0.1)
  expect_equal(2, sub1$theta, tolerance = 0.1)
  expect_equal(2, ipw1$theta, tolerance = 1)

  expect_equal(2, tml2$theta, tolerance = 0.1)
  expect_equal(2, sub2$theta, tolerance = 0.1)
  expect_equal(2, ipw2$theta, tolerance = 1)
})

test_that("For conditional estimation, trt_method must be riesz", {
  expect_error(
    lmtp_tmle(sim, trt = a, outcome = "Y", time_vary = time_vary,
              conditional = make_conditional(2), outcome_type = "continuous",
              shift = shift, folds = 1, riesz = FALSE)
  )
})
