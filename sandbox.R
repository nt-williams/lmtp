library(tidyverse)

n0 <- 1e6
data <- data.frame(matrix(NA, nrow = n0, ncol = 4))
names(data)  <- c(paste0('a', 1:2), paste0('y', 1:2))

for(t in 1:2) {
  if(t == 1) {
    data[, 'a1']  <- rbinom(n0, 1, 0.2)
    data[data$a1 == 0, 'y1'] <- rbinom(sum(data$a1 == 0), 1, 0.2)
    data[data$a1 == 1, 'y1'] <- 0
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
  }
}

mean(data[, 'y2'])

head(data)

at_risk_time_1 <- data$y1 == 0
m2 <- glm(y2 ~ a1 + a2, data = data, subset = at_risk_time_1)

q2_0 <- predict(m2, mutate(data, a2 = 0))
q2_0[!at_risk_time_1] <- 1

q2_1 <- predict(m2, mutate(data, a2 = 1))
q2_1[!at_risk_time_1] <- 1

m1_0 <- glm(q2_0 ~ a1, data = data)
m1_1 <- glm(q2_1 ~ a1, data = data)

1 - mean((1 - data$a1)*predict(m1_0, mutate(data, a1 = 0)) +
data$a1 * predict(m1_1, mutate(data, a1 = 0)))
