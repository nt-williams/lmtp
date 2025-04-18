# Add sporadic missingness to the first time point
sim_point_surv$Y.1[runif(nrow(sim_point_surv)) < 0.2 & sim_point_surv$C.1 == 1] <- NA_real_

A <- "trt"
Y <- paste0("Y.", 1:6)
C <- paste0("C.", 0:5)
W <- c("W1", "W2")

lmtp_curve(sim_point_surv, A, Y, W, cens = C, folds = 1,
           shift = static_binary_on, outcome_type = "survival")
