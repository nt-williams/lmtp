#
# set.seed(1234)
# n <- 200
# t_0 <- 6
# trt <- rbinom(n, 1, 0.5)
# adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n, 0, 2)))
# ftime <- round(1 + runif(n, 1, 4) - trt + adjustVars$W1 + adjustVars$W2)
# ftype <- round(runif(n, 0, 1))
# df <- data.frame(trt, adjustVars, ftime, ftype)
#
# wide_df <- prep_survival_data(Surv(ftime, ftype) ~ trt + W1 + W2, df, "trt", horizon = t_0)
#
# lmtp_tmle(wide_df$data, wide_df$trt, wide_df$outcome, wide_df$baseline,
#           cens = wide_df$cens, shift = static_binary_on)
#
# lmtp_tmle(wide_df$data, wide_df$trt, wide_df$outcome, wide_df$baseline,
#           cens = wide_df$cens, shift = static_binary_off)
#
# survtmle::survtmle(ftime = ftime, ftype = ftype,
#          trt = trt, adjustVars = adjustVars,
#          glm.ftime = "trt + W1 + W2",
#          method = "mean", t0 = t_0)
#
# temp <- data.frame(
#   id = c(1, 1, 2, 3, 4, 4, 5),
#   start = c(1, 3, 1, 1, 1, 4, 1),
#   stop = c(2, 4, 4, 2, 4, 5, 5),
#   baseline = c(3, 3, 1, 0, 2, 2, 1),
#   tv = c(0, 1, 1, 0, 1, 0, 1),
#   trt = c(2.5, 3, 4.6, 2.7, 1.2, 3.3, 4.4),
#   status = c(0, 0, 1, 0, 0, 1, 1)
# )
#
# prep_survival_data(Surv(start, stop, status) ~ trt + baseline + tv, temp,
#                    "trt", horizon = 5, id = "id")
#
# pf_temp <- unique(pf_temp)
#
# test <- prep_survival_data(Surv(start_day, stop_day, status) ~ pf_ratio + sex,
#                    pf_temp, "pf_ratio", horizon = 5, id = "patient_num")
