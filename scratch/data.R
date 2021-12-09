library(mvGPS)

n = 1e5

sim_dt <- gen_D(
  method = "u",
  n = n,
  rho_cond = 0.2,
  s_d1_cond = 2,
  s_d2_cond = 2,
  k = 3,
  C_mu = rep(0, 3),
  C_cov = 0.1,
  C_var = 1,
  d1_beta = c(0.5, 1, 0),
  d2_beta = c(0, 0.3, 0.75),
  seed = 06112020
)

D <- sim_dt$D
C <- sim_dt$C

alpha <- c(0.75, 1, 0.6, 1, 1)
sd_Y <- 1
X <- cbind(C, D)
Y <- X %*% alpha + rnorm(n, sd = sd_Y)

data = cbind(X, Y = Y)
data = as.data.frame(data)
names(data) = c(colnames(X), "Y")
multivariate_data = data

Dshift = D
Dshift[, 1] <- D[, 1] - 0.1
Dshift[, 2] <- D[, 2] - 0.5
Xshift <- cbind(C, Dshift)
Yshift <- Xshift %*% alpha + rnorm(n, sd = sd_Y)
