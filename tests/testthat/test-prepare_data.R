
set.seed(1234)
n <- 200
t_0 <- 6
trt <- rbinom(n, 1, 0.5)
adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
ftype <- round(runif(n, 0, 1))
df <- data.frame(trt, adjustVars, ftime, ftype)

survtmle(ftime = ftime, ftype = ftype,
         trt = trt, adjustVars = adjustVars,
         glm.ftime = "trt + W1 + W2",
         method = "mean", t0 = 6)

test <- prep_survival_data(df, "trt", "ftype", "ftime", covar = c("W1", "W2"), horizon = 6)

lmtp_surv_data(Surv(ftime, ftype) ~ trt + W1 + W2, df, "trt", horizon = 6)

lmtp_tmle(test$data, "trt", test$outcome_vars, c("W1", "W2"),
          cens = test$censoring_vars, shift = static_binary_on,
          outcome_type = "binomial", folds = 50)

lmtp_tmle(test$data, "trt", test$outcome_vars, c("W1", "W2"),
         cens = test$censoring_vars, shift = static_binary_off,
         outcome_type = "binomial", folds = 2)

