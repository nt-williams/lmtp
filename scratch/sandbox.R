library(tidyr)
library(dplyr)
library(mlr3superlearner)

head(sim_t4)

long =
  rename(sim_t4, Y_4 = Y) |>
  mutate(Y_3 = Y_4,
         Y_2 = Y_3,
         Y_1 = Y_2) |>
  pivot_longer(-ID, names_to = c(".value", "t"), names_sep = "_")

d = mutate(long, A = ifelse(A - 1 >= 1, A - 1, A))

long = mutate(long, across(c(L, A), as.factor))
d = mutate(d, across(c(L, A), as.factor))

f1 = mlr3superlearner(long[, c("t", "L", "A", "Y")],
                 "Y",
                 "ranger",
                 "binomial")

long$Y_f1 = predict(f1, d, type = "response")

mean(filter(long, t == 1)$Y_f1)

f2 = mlr3superlearner(long[long$t > 1, c("t", "L", "A", "Y_f1")],
                      "Y_f1",
                      "ranger",
                      outcome_type = "continuous")

long$Y_f2 = NA
long$Y_f2[long$t > 1] = predict(f2, filter(d, t > 1))

mean(filter(long, t == 2)$Y_f2)

f3 = mlr3superlearner(long[long$t > 2, c("t", "L", "A", "Y_f2")],
                      "Y_f2",
                      "ranger",
                      outcome_type = "continuous")

long$Y_f3 = NA
long$Y_f3[long$t > 2] = predict(f3, filter(d, t > 2))

mean(filter(long, t == 3)$Y_f3)

f4 = mlr3superlearner(long[long$t > 2, c("t", "L", "A", "Y_f3")],
                      "Y_f3",
                      "ranger",
                      outcome_type = "continuous")

long$Y_f4 = NA
long$Y_f4[long$t > 3] = predict(f4, filter(d, t > 3))

mean(filter(long, t == 4)$Y_f4)
