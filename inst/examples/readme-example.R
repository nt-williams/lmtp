
library(lmtp)

# A data set with treatment nodes at 4 time points and a binary outcome at time 5
head(sim_t4)

# Define treatment variables, covariates, and a treatment policy
a <- c("A_1", "A_2", "A_3", "A_4")
nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

d <- function(A) {
  delta <- 1
  return((A - delta) * (A - delta >= 0) + A * (A - delta < 0))
}

# Define a set of sl3 learners
lrnrs <- sl3::make_learner_stack(sl3::Lrnr_glm, sl3::Lrnr_mean)

# calculate estimates!
lmtp_tmle(sim_t4, a, "Y", nodes, k = 0, shift = d, outcome_type = "binomial",
          learners_outcome = lrnrs, learners_trt = lrnrs)
