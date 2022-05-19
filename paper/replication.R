# Replication script for the lmtp package Observational Studies paper

# Nicholas Williams, MPH

# the following code may be used to install lmtp and other necessary packages
# install.packages(c("lmtp", "devtools", "earth", "ranger"))

library(lmtp)

# Section: 3.2 Required data structure ------------------------------------

W <- c("W_1", "W_2")
A <- c("A_1", "A_2")
L <- list(c("L_11", "L_12"), c("L_21", "L_22"))
create_node_list(trt = A, baseline = W, time_vary = L, tau = 2)

# Example 1 ---------------------------------------------------------------

A <- c("A_1", "A_2", "A_3", "A_4")
L <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

f <- function(data, trt) {
    (data[[trt]] - 1) * (data[[trt]] - 1 >= 1) +
        data[[trt]] * (data[[trt]] - 1 < 1)
}

set.seed(4524)
lmtp_tmle(sim_t4, A, "Y", time_vary = L, shift = f,
          intervention_type = "mtp", folds = 1)

# Example 2 ---------------------------------------------------------------

A <- c("A1", "A2")
C <- c("C1", "C2")
L <- list(c("L1"), c("L2"))

f <- function(data, trt) data[[trt]] + 0.5
sl_lib <- c("SL.glm", "SL.ranger", "SL.earth")

set.seed(624)
lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = C,
         shift = f, intervention_type = "mtp", folds = 1,
         learners_trt = sl_lib, learners_outcome = sl_lib)

lmtp_sdr(sim_cens, A, "Y", time_vary = L, cens = C,
         shift = NULL, folds = 1, learners_trt = sl_lib, learners_outcome = sl_lib)

# Example 3 ---------------------------------------------------------------

A <- "trt"
W <- c("W1", "W2")
C <- paste0("C.", 0:5)
Y <- paste0("Y.", 1:6)

set.seed(5423)
a_1 <- lmtp_tmle(sim_point_surv, A, Y, W, cens = C,
                 outcome_type = "survival", folds = 1,
                 learners_trt = sl_lib,
                 learners_outcome = sl_lib,
                 shift = static_binary_on)

a_0 <- lmtp_tmle(sim_point_surv, A, Y, W, cens = C,
                 outcome_type = "survival", folds = 1,
                 learners_trt = sl_lib,
                 learners_outcome = sl_lib,
                 shift = static_binary_off)

lmtp_contrast(a_1, ref = a_0, type = "rr")

# Example 4 ---------------------------------------------------------------

A <- c("A_1", "A_2", "A_3", "A_4")
L <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

f <- function(data, trt) {
  mtp <- function(data, trt) {
    (data[[trt]] - 1) * (data[[trt]] - 1 >= 1) +
      data[[trt]] * (data[[trt]] - 1 < 1)
  }

  if (trt == "A_1") {
    return(mtp(data, trt))
  }

  ifelse(
    data[[sub("A", "L", trt)]] == 1,
    mtp(data, trt),
    data[[trt]]
  )
}

set.seed(697)
lmtp_sdr(sim_t4, A, "Y", time_vary = L, shift = f,
         intervention_type = "mtp", folds = 1,
         learners_outcome = sl_lib, learners_trt = sl_lib)
