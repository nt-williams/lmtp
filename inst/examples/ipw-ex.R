\dontrun{

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
  psi1.1 <- lmtp_ipw(ex1_dat, "A", "Y", list(c("W")), shift = d, folds = 2)
  psi1.1

  # Example 1.2
  # Point treatment, continuous exposure, continuous outcome, no loss-to-follow-up
  # Interesed in the effect of a modified treatment policy where A is decreased by 15
  # units only among observations whose observed A was above 80.
  # The true value under this intervention is about 513.
  d <- function(x) (x > 80)*(x - 15) + (x <= 80)*x
  psi1.2 <- lmtp_ipw(ex1_dat, "A", "Y", list(c("W")), shift = d, folds = 2)
  psi1.2

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
    psi2.1 <- lmtp_ipw(sim_t4, a, "Y", time_varying, shift = d, folds = 2)
  })
  psi2.1

  # Example 2.2
  # Example 2.1 assumed that the outcome (as well as the treatment variables)
  # were directly affected by all other nodes in the past. In certain situtations,
  # domain specific knowledge may suggest otherwise leading to a Markov processes.
  # This can be controlled using the k argument.
  progressr::with_progress({
    psi2.2 <- lmtp_ipw(sim_t4, a, "Y", time_varying, shift = d,
                       k = 1, folds = 2)
  })
  psi2.2

  # Example 2.3
  # Using the same data as examples 2.1 and 2.3, but now treating the exposure
  # as an ordered categorical variable. To account for the exposure being a
  # factor we just need to modify the shift function (and the original data)
  # so as to respect this.
  for (i in a) {
    sim_t4[[i]] <- factor(sim_t4[[i]], levels = 0:5, ordered = T)
  }

  d <- function(a) {
    out <- list()
    for (i in 1:length(a)) {
      if (as.character(a[i]) %in% c("0", "1")) {
        out[[i]] <- as.character(a[i])
      } else {
        out[[i]] <- as.numeric(as.character(a[i])) - 1
      }
    }
    factor(unlist(out), levels = 0:5, ordered = T)
  }

  progressr::with_progress({
    psi2.3 <- lmtp_ipw(sim_t4, a, "Y", time_varying, shift = d, folds = 2)
  })
  psi2.3

  # Example 3.1
  # Longitudinal setting, time-varying binary treatment, time-varying covariates
  # and baseline covariates with no loss-to-follow-up. Interested in a traditional
  # causal effect where treatment is set to 1 at all time points for all observations.
  data("iptwExWide", package = "twang")
  a <- paste0("tx", 1:3)
  baseline <- c("gender", "age")
  nodes <- list(c("use0"), c("use1"), c("use2"))
  progressr::with_progress({
    psi3.1 <-
      lmtp_ipw(iptwExWide, a, "outcome", nodes, baseline = baseline,
               shift = function(x) 1, folds = 2)
  })
  psi3.1

  # Example 4.1
  # Longitudinal setting, time-varying continuous treatment, time-varying covariates,
  # binary outcome with right censoring. Interested in the mean population outcome under
  # the observed exposures in a hypothetical population with no loss-to-follow-up.
  head(sim_cens)
  a <- c("A1", "A2")
  nodes <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  y <- "Y"
  psi4.1 <- lmtp_ipw(sim_cens, a, y, nodes, cens = cens, shift = NULL, folds = 2)
  psi4.1

  # Example 4.2
  # Using the same data as example 4.1, but now interested in the causal effect of a
  # treatment policy where exposure increased by 0.5 units at all time points. The
  # true value under this intervention is about 0.88.
  d <- function(x) x + 0.5
  psi4.2 <- lmtp_ipw(sim_cens, a, y, nodes, cens = cens, shift = d, folds = 2)
  psi4.2

  # Example 5.1
  # Time-to-event analysis with a binary time-invariant exposure. Interested in
  # the effect of treatment being given to all observations on the cumulative
  # incidence of our time-to-event outcome.
  a <- "trt"
  # for a survival problem, the outcome arugment now takes a vector of outcomes
  # if an observation experiences the event prior to the end of follow-up, all future
  # outcome nodes should be set to 1 (i.e., last observation carried forward).
  y <- paste0("Y", 0:6)
  cens <- paste0("C.", 0:5)
  baseline <- c("W1", "W2")
  # even if there are no time varying covariates, we still provide a list the same
  # length as the number of time points, but with items specified as NULL.
  nodes <- lapply(0:5, function(x) NULL)
  progressr::with_progress({
    psi5.1 <- lmtp_ipw(sim_point_surv, a, y, nodes, baseline, cens,
                       shift = function(x) 1, folds = 2)
  })
  psi5.1

}
