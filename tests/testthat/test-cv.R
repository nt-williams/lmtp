
library(sl3)

a <- c("A_1", "A_2", "A_3", "A_4")
nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

d <- function(a) (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)

glm_learner <- Lrnr_glm$new()
mean_learner <- Lrnr_mean$new()
lrnrs <- Stack$new(glm_learner, mean_learner)

lmtp_ipw(sim_t4, a, "Y", nodes, k = 1, shift = d,
         learners = lrnrs, V = 5, progress_bar = F)

lmtp_sub(sim_t4, a, "Y", nodes, k = 1, shift = d,
         learners = lrnrs, V = 5, progress_bar = F)

lmtp_tmle(sim_t4, a, "Y", nodes, k = 1, shift = d,
         learners_outcome = lrnrs, learners_trt = lrnrs,
         V = 5, progress_bar = F)


meta <- Meta$new(
  sim_t4,
  a,
  "Y",
  nodes,
  NULL,
  NULL,
  k = 1,
  shift = d,
  V = 4
)

estimate_tmle(meta$data[[4]]$train, meta$shifted_data[[4]]$train,
              meta$data[[4]]$valid, meta$shifted_data[[4]]$valid,
              "xyz", meta$node_list, NULL, meta$tau, meta$tau,
              "binomial", meta$m[[3]], meta$m[[3]], dr,
              make_learner(Lrnr_glm), check_pb(FALSE, meta$tau, "Estimating propensity"))

estimate_sub(
  meta$data[[1]]$train,
  meta$shifted_data[[1]]$train,
  meta$shifted_data[[1]]$valid,
  "Y",
  meta$node_list,
  NULL,
  meta$tau,
  "binomial",
  sl3::make_learner(sl3::Lrnr_glm),
  meta$m[[1]]$valid,
  pb = check_pb(FALSE, meta$tau, ".")
)

ce <-
  estimate_c(sim_t4, meta$data[[4]]$train, meta$data[[4]]$valid, NULL, "Y", meta$tau,
             meta$node_list, sl3::make_learner(sl3::Lrnr_glm))

dr <- estimate_r(
  training = meta$data[[4]]$train,
  validation = meta$data[[4]]$valid,
  trt = c("A_1", "A_2", "A_3", "A_4"),
  cens = NULL,
  C = ce,
  shift = d,
  tau = meta$tau,
  node_list = meta$node_list,
  learners = sl3::make_learner(sl3::Lrnr_glm),
  pb = check_pb(FALSE, meta$tau, "Estimating propensity")
) %>% ratio_ite()
