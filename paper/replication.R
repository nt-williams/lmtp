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

# Example 3 ---------------------------------------------------------------

x <- "trt"
w <- c("W1", "W2")
cen <- paste0("C.", 0:5)
y <- paste0("Y.", 1:6)

set.seed(5423)
tml1 <- lmtp_tmle(sim_point_surv, x, y, w, cens = cen,
                  learners_trt = lrnrs, learners_outcome = lrnrs,
                  shift = static_binary_on, folds = 3)

set.seed(6354)
tml0 <- lmtp_tmle(sim_point_surv, x, y, w, cens = cen,
                  learners_trt = lrnrs, learners_outcome = lrnrs,
                  shift = static_binary_off, folds = 3)

lmtp_contrast(tml1, ref = tml0, type = "rr")

set.seed(432)
sdr1 <- lmtp_sdr(sim_point_surv, x, y, w, cens = cen,
                 learners_trt = lrnrs, learners_outcome = lrnrs,
                 shift = static_binary_on, folds = 3)

set.seed(143)
sdr0 <- lmtp_sdr(sim_point_surv, x, y, w, cens = cen,
                 learners_trt = lrnrs, learners_outcome = lrnrs,
                 shift = static_binary_off, folds = 3)

lmtp_contrast(sdr1, ref = sdr0, type = "rr")

# Example 4 ---------------------------------------------------------------

x <- c("A_1", "A_2", "A_3", "A_4")
tv <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
y <- "Y"

shift <- function(data, trt) {
  (data[[trt]] - 1) * (data[[trt]] - 1 >= 1) +
    data[[trt]] * (data[[trt]] - 1 < 1)
}

dynamic_shift <- function(data, trt) {
  if (trt == "A_1") {
    shift(data, trt)
  } else {
    ifelse(data[[sub("A", "L", trt)]] == 1,
           shift(data, trt),
           data[[trt]])
  }
}

set.seed(540)
lmtp_tmle(sim_t4, x, y, time_vary = tv, shift = dynamic_shift, folds = 3,
          learners_outcome = lrnrs, learners_trt = lrnrs)

set.seed(697)
lmtp_sdr(sim_t4, x, y, time_vary = tv, shift = dynamic_shift, folds = 3,
         learners_outcome = lrnrs, learners_trt = lrnrs)
