library(lmtp)
# Example 1.1
# Point treatment, continuous exposure, continuous outcome, no loss-to-follow-up
# Interested in the effect of a population wide decrease in A of 5 units
# The true value under this intervention is about 519.
set.seed(56)
n <- 1000
W <- rnorm(n, 10, 5)
A <- 23 + 5*W + rnorm(n)
Y <- 7.2*A + 3*W + rnorm(n)
ex1_dat <- data.frame(W, A, Y)
d <- function(x) x - 5
fit1.1 <- lmtp_tmle(ex1_dat, "A", "Y", list(c("W")), shift = d,
                    outcome_type = "continuous", folds = 2)

# Example 1.2
# Point treatment, continuous exposure, continuous outcome, no loss-to-follow-up
# Interesed in the effect of a modified treatment policy where A is decreased by 15
# units only among observations whose observed A was above 80.
# The true value under this intervention is about 513.
d <- function(x) (x > 80)*(x - 15) + (x <= 80)*x
fit1.2 <- lmtp_tmle(ex1_dat, "A", "Y", list(c("W")), shift = d,
                    outcome_type = "continuous", folds = 2)
fit1.2

# Example 2.1
# Longitudinal setting, time-varying continuous exposure bounded by 0,
# time-varying covariates, and a binary outcome with no loss-to-follow-up.
# Interested in the effect of a treatment policy where exposure decreases by
# one unit at every time point if an observations observed exposure is greater
# than or equal to 2. The true value under this intervention is about 0.305.
head(sim_t4)
# specifying treament variables
a <- c("A_1", "A_2", "A_3", "A_4")
# specifying time varying covariates
time_varying <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
# treatment policy function to be applied at all time points
d <- function(a) {
  (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)
}
progressr::with_progress({
  fit2.1 <- lmtp_tmle(sim_t4, a, "Y", time_varying, shift = d, folds = 2)
})
fit2.1

# Example 2.2
# Example 2.1 assumed that the outcome (as well as the treatment variables)
# were directly affected by all other nodes in the past. In certain situtations,
# domain specific knowledge may suggest otherwise leading to a Markov processes.
# This can be controlled using the k argument.
progressr::with_progress({
  fit2.2 <- lmtp_tmle(sim_t4, a, "Y", time_varying, shift = d,
                      k = 1, folds = 2)
})
fit2.2

# Example 3.1
# Longitudinal setting, time-varying binary treatment, time-varying covariates
# and baseline covariates with no loss-to-follow-up. Interested in a traditional
# causal effect where treatment is set to 1 at all time points for all observations.
data("iptwExWide", package = "twang")
a <- paste0("tx", 1:3)
baseline <- c("gender", "age")
nodes <- list(c("use0"), c("use1"), c("use2"))
progressr::with_progress({
  fit3.1 <-
    lmtp_tmle(iptwExWide, a, "outcome", nodes, baseline = baseline,
              shift = function(x) 1, outcome_type = "continuous",
              folds = 2)
})
fit3.1



