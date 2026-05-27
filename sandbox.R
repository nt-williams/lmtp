library(tidyverse)

set.seed(435634)

n0 <- 1e3
data <- data.frame(matrix(NA, nrow = n0, ncol = 9))
names(data)  <- c(paste0('a', 1:3), paste0('y', 1:3), paste0("c", 1:3))

for(t in 1:3) {
  if(t == 1) {
    data[, 'a1']  <- rbinom(n0, 1, 0.2)
    data[data$a1 == 0, 'y1'] <- rbinom(sum(data$a1 == 0), 1, 0.2)
    data[data$a1 == 1, 'y1'] <- 0
    data[, 'c1'] <- rbinom(n0, 1, 0.9)  # MCAR: c=1 still observed, c=0 lost to follow-up
  } else {
    vaccinated  <- rowSums(data[, paste0('a', 1:(t-1)), drop = FALSE]) > 0
    alive  <- rowSums(data[, paste0('y', 1:(t-1)), drop = FALSE]) == 0
    data[!vaccinated & alive, paste0('a', t)] <- rbinom(sum(!vaccinated & alive), 1, 0.2)
    data[ vaccinated & alive, paste0('a', t)] <- 1
    data[!alive, paste0('a', t)]  <- vaccinated[!alive]

    vaccinated  <- rowSums(data[, paste0('a', 1:t), drop = FALSE]) > 0
    data[ vaccinated & alive, paste0('y',t)] <- 0
    data[!vaccinated & alive, paste0('y',t)] <- rbinom(sum(!vaccinated & alive), 1, 0.2)
    data[!alive, paste0('y',t)] <- 1

    still_observed <- data[, paste0("c", t-1)] == 1
    data[ still_observed, paste0("c", t)] <- rbinom(sum(still_observed), 1, 0.9)  # MCAR
    data[!still_observed, paste0("c", t)] <- 0
  }
}

data$y2[data$c1 == 0] <- NA_real_
data$y1[data$c1 == 0] <- NA_real_
data$a2[data$c1 == 0] <- NA_real_
data$y3[data$c1 == 0] <- NA_real_
data$a3[data$c1 == 0] <- NA_real_

data$y2[data$c2 == 0] <- NA_real_
data$y3[data$c2 == 0] <- NA_real_
data$a3[data$c2 == 0] <- NA_real_

data$y3[data$c3 == 0] <- NA_real_

# 1 - mean(data[, 'y2'])
#
# head(data)
#
# at_risk_time_1 <- data$y1 == 0
# m2 <- glm(y2 ~ a1 + a2, data = data, subset = at_risk_time_1)
#
# augmented <- delay_augment.data.frame(data, task$sequences(1))
# augmented$a2 <- one_time_delay(augmented, c("..i..lmtp_tmp_s1", "a2"))
#
# q2 <- predict(m2, augmented)
# q2[augmented$y1 == 1] <- 1
#
# m1 <- glm(q2 ~ ..i..lmtp_tmp_s1*a1, data = augmented)
#
# q1 <- predict(m1, mutate(data, ..i..lmtp_tmp_s1 = a1, a1 = 0))
#
# 1 - mean(q1)

# testing functions -------------------------------------------------------

# data$w <- rnorm(nrow(data))
data$w <- rep(1, nrow(data))
trt <- c("a1", "a2", "a3")
baseline <- list(
  trt     = "w",
  cens    = "w",
  outcome = "w"
)
time_vary <- NULL
cens <- c("c1", "c2", "c3")
compete <- NULL
id <- NULL
outcome <- c("y1", "y2", "y3")

# variable_names <- c(unlist(trt), outcome, unlist(time_vary), baseline, cens, compete, id)

# task <- LmtpTask$new(
#   data = data,
#   shifted = NULL,
#   A = trt,
#   Y = outcome,
#   L = time_vary,
#   W = baseline,
#   C = cens,
#   D = compete,
#   k = Inf, id = id,
#   outcome_type = "survival",
#   folds = 10,
#   weights = NULL,
#   bounds = NULL
# )

# prop <- cf_propensity_score(task, "SL.glm", lmtp_control(), NULL)

# cf_tmle_delay(task, prop$propensity_score, "SL.glm.interaction", lmtp_control(), NULL)

# foo <- estimate_sdr_delay(task, 1, prop$propensity_score, "SL.glm.interaction", lmtp_control(), NULL)
# ife::ife(mean(foo$predictions[[1]]), foo$efficient_influence_function)


# htlmtp_sdr(data, trt, outcome, baseline,
#            outcome_type = "survival",
#            learners_outcome = "SL.glm.interaction",
#            learners_trt = "SL.glm",
#            learners_cens = "SL.glm",
#            folds = 10,
#            control = lmtp_control())


wc <- htlmtp_sdr(data, trt, outcome, baseline,
            cens = cens,
           outcome_type = "survival",
           learners_outcome = "SL.glm.interaction",
           learners_trt = "SL.glm",
           learners_cens = "SL.glm",
           folds = 1,
           control = lmtp_control())

noc <- htlmtp_tmle(data, trt, outcome, baseline,
           outcome_type = "survival",
           learners_outcome = "SL.glm.interaction",
           learners_trt = "SL.glm",
           learners_cens = "SL.glm",
           folds = 1,
           control = lmtp_control())
