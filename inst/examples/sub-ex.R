\donttest{
  set.seed(56)
  n <- 1000
  W <- rnorm(n, 10, 5)
  A <- 23 + 5*W + rnorm(n)
  Y <- 7.2*A + 3*W + rnorm(n)
  ex1_dat <- data.frame(W, A, Y)

  # Example 1.1
  # Point treatment, continuous exposure, continuous outcome, no loss-to-follow-up
  # Interested in the effect of a modified treatment policy where A is decreased by 15
  # units only among observations whose observed A was above 80.
  # The true value under this intervention is about 513.
  d <- function(data, x) (data[[x]] > 80)*(data[[x]] - 15) + (data[[x]] <= 80)*data[[x]]
  lmtp_sub(ex1_dat, "A", "Y", "W", shift = d, outcome_type = "continuous", folds = 2)

  # Example 2.1
  # Longitudinal setting, time-varying continuous exposure bounded by 0,
  # time-varying covariates, and a binary outcome with no loss-to-follow-up.
  # Interested in the effect of a treatment policy where exposure decreases by
  # one unit at every time point if an observations observed exposure is greater
  # than or equal to 2. The true value under this intervention is about 0.305.
  head(sim_t4)
  a <- c("A_1", "A_2", "A_3", "A_4")
  tv <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
  d <- function(data, trt) {
    a <- data[[trt]]
    (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)
  }

  # BONUS: progressr progress bars!
  progressr::handlers(global = TRUE)

  lmtp_sub(sim_t4, a, "Y", time_vary = tv, shift = d, folds = 2)

  # Example 2.2
  # The previous example assumed that the outcome (as well as the treatment variables)
  # were directly affected by all other nodes in the past. In certain situations,
  # domain specific knowledge may suggest otherwise.
  # This can be controlled using the k argument.
  lmtp_sub(sim_t4, a, "Y", time_vary = tv, shift = d, k = 0, folds = 2)

  # Example 2.3
  # Using the same data as examples 2.1 and 2.2.
  # Now estimating the effect of a dynamic modified treatment policy.
  a <- c("A_1", "A_2", "A_3", "A_4")
  time_varying <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

  # creating a dynamic mtp that applies the shift function
  # but also depends on history and the current time
  dynamic_mtp <- function(data, trt) {
    mtp <- function(data, trt) {
      (data[[trt]] - 1) * (data[[trt]] - 1 >= 1) + data[[trt]] * (data[[trt]] - 1 < 1)
    }

    # if its the first time point, follow the same mtp as before
    if (trt == "A_1") return(mtp(data, trt))

    # otherwise check if the time varying covariate equals 1
    ifelse(
      data[[sub("A", "L", trt)]] == 1,
      mtp(data, trt), # if yes continue with the policy
      data[[trt]]     # otherwise do nothing
    )
  }
  psi2.3 <- lmtp_sub(sim_t4, a, "Y", time_vary = time_varying,
                     k = 0, shift = dynamic_mtp, folds = 2)
  psi2.3

  # Example 2.4
  # Using the same data as examples 2.1, 2.2, and 2.3, but now treating the exposure
  # as an ordered categorical variable. To account for the exposure being a
  # factor we just need to modify the shift function (and the original data)
  # so as to respect this.
  tmp <- sim_t4
  for (i in a) {
    tmp[[i]] <- factor(tmp[[i]], levels = 0:5, ordered = TRUE)
  }
  d <- function(data, trt) {
    out <- list()
    a <- data[[trt]]
    for (i in 1:length(a)) {
      if (as.character(a[i]) %in% c("0", "1")) {
        out[[i]] <- as.character(a[i])
      } else {
        out[[i]] <- as.numeric(as.character(a[i])) - 1
      }
    }
    factor(unlist(out), levels = 0:5, ordered = TRUE)
  }
  lmtp_sub(tmp, a, "Y", time_vary = tv, shift = d, k = 0, folds = 2)

  # Example 3.1
  # Longitudinal setting, time-varying binary treatment, time-varying covariates
  # and baseline covariates with no loss-to-follow-up. Interested in a "traditional"
  # causal effect where treatment is set to 1 at all time points for all observations.
  if (require("twang")) {
    data("iptwExWide", package = "twang")
    a <- paste0("tx", 1:3)
    baseline <- c("gender", "age")
    tv <- list(c("use0"), c("use1"), c("use2"))
    lmtp_sub(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
              shift = static_binary_on, outcome_type = "continuous")
  }

  # Example 4.1
  # Longitudinal setting, time-varying continuous treatment, time-varying covariates,
  # binary outcome with right censoring. Interested in the mean population outcome under
  # the observed exposures in a hypothetical population with no loss-to-follow-up.
  head(sim_cens)
  a <- c("A1", "A2")
  tv <- list(c("L1"), c("L2"))
  cens <- c("C1", "C2")
  y <- "Y"
  lmtp_sub(sim_cens, a, y, time_vary = tv, cens = cens, shift = NULL, folds = 2)

  # Example 5.1
  # Time-to-event analysis with a binary time-invariant exposure. Interested in
  # the effect of treatment being given to all observations on the cumulative
  # incidence of the outcome.
  # For a survival problem, the outcome argument now takes a vector of outcomes
  # if an observation experiences the event prior to the end of follow-up, all future
  # outcome nodes should be set to 1 (i.e., last observation carried forward).
  a <- "trt"
  y <- paste0("Y.", 1:6)
  cens <- paste0("C.", 0:5)
  baseline <- c("W1", "W2")
  lmtp_sub(sim_point_surv, a, y, baseline, cens = cens, folds = 2,
            shift = static_binary_on, outcome_type = "survival")
}
