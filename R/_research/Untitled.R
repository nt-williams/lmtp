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

# shifted1 <- data.frame(a1 = rep(0, n0), a2 = dat$a1)
shiftedDelayed <- data.frame(a1 = rep(0, n0), a2 = data$a1)
shiftedStochastic <- data.frame(a1 = rep(0, n0), a2 = as.numeric(runif(n0) < 0.2))

fshift <- function(shifted){
  for(t in 1:2) {
    if(t == 1) {
      shifted[shifted$a1 == 0, 'y1'] <- rbinom(sum(shifted$a1 == 0), 1, 0.2)
      shifted[shifted$a1 == 1, 'y1'] <- 0
    } else {
      alive  <- rowSums(shifted[, paste0('y', 1:(t-1)), drop = FALSE]) == 0
      vaccinated  <- rowSums(shifted[, paste0('a', 1:t), drop = FALSE]) > 0
      shifted[ vaccinated & alive, paste0('y',t)] <- 0
      shifted[!vaccinated & alive, paste0('y',t)] <- rbinom(sum(!vaccinated & alive), 1, 0.2)
      shifted[!alive, paste0('y',t)] <- 1
    }
  }
  shifted <- shifted[, c(paste0('a', 1:2), paste0('y', 1:2))]
  return(shifted)
}

shiftedDelayed <- fshift(shiftedDelayed)
shiftedStochastic <- fshift(shiftedStochastic)

## True values
1 - mean(shiftedDelayed[, 'y2'])
1 - mean(shiftedStochastic[, 'y2'])

# G-formula
m3_d <- data$y2

at_risk_time_1 <- data$y1 == 0
m2 <- glm(m3_d ~ a1 + a2, data = data, subset = at_risk_time_1)
data$m2_d <- NA_real_
data$m2_d[at_risk_time_1] <- predict(m2, mutate(data, a2 = shiftedDelayed$a2)[at_risk_time_1, ])
data$m2_d[!at_risk_time_1] <- 1

m1 <- glm(m2_d ~ a1, data = data)
data$m1_d <- predict(m1, mutate(data, a1 = shiftedDelayed$a1))

1 - mean(data$m1_d)
