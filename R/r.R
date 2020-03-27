
# set.seed(429153)
# n_obs <- 50
# W <- replicate(2, rbinom(n_obs, 1, 0.5))
# A <- rnorm(n_obs, mean = 2 * W, sd = 1)
# Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
# tl <- list(c("X1", "X2", "A"))
# df <- data.frame(W, A, Y)
# us <- shift_data(df, "A", 0.5)
# ds <- shift_data(df, "A", -0.5)
# r <- matrix(nrow = n_obs, ncol = 1)

estimate_r <- function(data, A, shifted_up, shifted_down,
                       tau = NULL, node_list = NULL, r = NULL) {

  n <- nrow(data)
  a <- data[, A]
  us <- shifted_up[, A]
  ds <- shifted_down[, A]
  w <- data[, c("X1", "X2")]
  gn <- haldensify::haldensify(A = a, W = w, n_bins = c(5, 10))

  pred <- data.frame(obs = rep(NA, n), upshift = rep(NA, n), downshift = rep(NA, n))
  pred$obs <- stats::predict(object = gn, new_A = a, new_W = w)
  pred$upshift <- stats::predict(object = gn, new_A = us, new_W = w)
  pred$downshift <- stats::predict(object = gn, new_A = ds, new_W = w)

  # TODO compute some ratios
}
