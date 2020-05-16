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
progressr::with_progress({
  fit1.1 <- lmtp_tmle(ex1_dat, "A", "Y", list(c("W")), shift = d,
                      outcome_type = "continuous", folds = 2)
})

# Example 1.2
# Point treatment, continuous exposure, continuous outcome, no loss-to-follow-up
# Interesed in the effect of a modified treatment policy where A is decreased by 15
# units only among observations whose observed A was above 80.
# The true value under this intervention is about 513.
d <- function(x) (x > 80)*(x - 15) + (x <= 80)*x
progressr::with_progress({
  fit1.2 <- lmtp_tmle(ex1_dat, "A", "Y", list(c("W")), shift = d,
                      outcome_type = "continuous", folds = 2)
})

