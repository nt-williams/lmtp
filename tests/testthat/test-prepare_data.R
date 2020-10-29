
set.seed(1234)
n <- 200
t_0 <- 6
trt <- rbinom(n, 1, 0.5)
adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
ftype <- round(runif(n, 0, 1))
test <- data.frame(trt, adjustVars, ftime, ftype)

test <- prepare_data(test, "trt", "ftype", "ftime")

lmtp_tmle(test$data, "trt", test$outcome_vars, c("W1", "W2"),
          cens = test$censoring_vars, shift = static_binary_on,
          outcome_type = "binomial")

library(survrct)

surv <- survrct(Surv(ftime, ftype) ~ trt + W1 + W2,
                target = "trt", data = test, estimator = "tmle")
